"""
Code to generate bps yaml files for single frame processing.
"""
import os
import sqlite3
import numpy as np
import pandas as pd

__all__ = ['SfpYamlFactory']

# Template for making bps yaml files.  Use '$(...)' for env vars to be
# resolved by bps.  The '$()' will be replaced by '${}' after the
# local python variables have been inserted.
BPS_SFP_YAML = """includeConfigs:
  - $(GEN3_WORKFLOW_DIR)/python/desc/gen3_workflow/etc/bps_drp_baseline.yaml
  - $(GEN3_WORKFLOW_DIR)/examples/bps_DC2-3828-y1_resources.yaml

pipelineYaml: "$(OBS_LSST_DIR)/pipelines/imsim/DRP.yaml#singleFrame"

payload:
  inCollection: LSSTCam-imSim/defaults
  payloadName: {payloadName}
  butlerConfig: {repo}
  dataQuery: "{dataQuery}"
"""

class SfpYamlFactory:
    """
    Class to generate bps yaml files for single frame processing given
    an sqlite3 db of the overlaps between ccd-visits and skymap tracts
    and patches.
    """
    def __init__(self, overlap_db, repo):
        """
        Parameters
        ----------
        overlap_db : sqlite3 db file
            This file contains an 'overlaps' table with the ccd-visits
            that overlap each patch in the repo skymap.
        """
        if not os.path.isfile(overlap_db):
            raise FileNotFoundError(f'{overlap_db} not found')
        self.overlap_db = overlap_db
        self.repo = repo

    def create(self, tracts, num_parts=1, visit_range=None,
               processed_visits=None):
        """
        Create bps yaml files for single-frame processing of the
        visits that overlap with the specified tracts.

        Parameters
        ----------
        tracts : list-like
            List of tract numbers in the repo skymap to consider.
        num_parts : int [1]
            Number of parts to divide processing among the corresponding
            number of bps yaml files.
        visit_range : (int, int) [None]
            Range of visits to consider.
        processed_visits : list [None]
            List of visits that have already been processed and which
            should be excluded.

        Returns
        -------
        list of visits included in the bps yaml files.
        """
        repo = self.repo
        if num_parts < 1:
            raise ValueError('Must have num_parts >= 1.')
        tract_list = ','.join([str(_) for _ in tracts])
        query = f'select * from overlaps where tract in ({tract_list})'
        if visit_range is not None:
            query += (f' and visit >= {visit_range[0]}'
                      f' and visit <= {visit_range[1]}')
        if processed_visits is not None and processed_visits:
            processed_visit_list \
                = '(' + ','.join([str(_) for _ in processed_visits]) + ')'
            query += f' and visit not in {processed_visit_list}'
        with sqlite3.connect(self.overlap_db) as con:
            df0 = pd.read_sql(query, con)
        visits = sorted(list(set(df0['visit'])))

        indexes = np.linspace(0, len(visits) + 1, num_parts + 1, dtype=int)

        tract_list = '_'.join([str(_) for _ in tracts])
        for part, (imin, imax) in enumerate(zip(indexes[:-1], indexes[1:])):
            # Use '[...]' here, and replace with '(...)' after env var
            # delimiters have been replaced.
            visit_list = '[' + ','.join([str(_) for _ in visits[imin:imax]]) \
                         + ']'
            payloadName = f'sfp_Y1_{tract_list}_visits_part_{part:02d}'
            dataQuery = f"instrument='LSSTCam-imSim' and visit in {visit_list}"
            outfile = f'bps_{payloadName}.yaml'
            print(outfile)
            with open(outfile, 'w') as output:
                output_string \
                    = BPS_SFP_YAML.format(**locals())\
                                  .replace('(', '{').replace(')', '}')\
                                  .replace('[', '(').replace(']', ')')
                output.write(output_string)

        return visits
