import os
import sqlite3
import numpy as np
import pandas as pd

overlap_db = ('/global/cfs/cdirs/lsst/shared/DC2-prod/Run2.2i/'
              'desc_dm_drp/v19.0.0/rerun/run2.2i-calexp-v1/'
              'tracts_mapping.sqlite3')

#dia_tracts =  (4430, 4431, 4432, 4638, 4639, 4640)
tract = 4430

assert os.path.isfile(overlap_db)

with sqlite3.connect(overlap_db) as con:
    df0 = pd.read_sql(f'select * from overlaps where tract={tract}', con)

y1_visit_max = 262622

df = df0.query(f'visit <= {y1_visit_max}')

visits = sorted(list(set(df['visit'])))
print('total number of visits:', len(visits))
indexes = np.linspace(0, len(visits) + 1, 2, dtype=int)
print(indexes)

repo = '/global/cfs/cdirs/lsst/production/gen3/DC2/Run2.2i/repo'

nodes = 10
memory = 90000

# template for making bps yaml files
bps_sfp_yaml = """includeConfigs:
  - $(GEN3_WORKFLOW_DIR)/python/desc/gen3_workflow/etc/bps_drp_baseline.yaml
  - $(GEN3_WORKFLOW_DIR)/examples/bps_DC2-3828-y1_resources.yaml

pipelineYaml: "$(OBS_LSST_DIR)/pipelines/imsim/DRP.yaml#singleFrame"

payload:
  inCollection: LSSTCam-imSim/defaults
  payloadName: {payloadName}
  butlerConfig: {repo}
  dataQuery: "{dataQuery}"

parsl_config:
  retries: 1
  monitoring: true
  executor: WorkQueue
  provider: Local
  nodes_per_block: {nodes}
  worker_options: "--memory={memory}"
"""

for part, (imin, imax) in enumerate(zip(indexes[:-1], indexes[1:])):
    visit_list = '[' + ','.join([str(_) for _ in visits[imin:imax]]) + ']'
    payloadName = f'sfp_Y1_{tract}_visits_part_{part:02d}'
    dataQuery = f"instrument='LSSTCam-imSim' and visit in {visit_list}"
    outfile = f'bps_{payloadName}.yaml'
    print(outfile)
    with open(outfile, 'w') as output:
        output_string = bps_sfp_yaml.format(**locals())\
                                    .replace('(', '{').replace(')', '}')\
                                    .replace('[', '(').replace(']', ')')
        output.write(output_string)
