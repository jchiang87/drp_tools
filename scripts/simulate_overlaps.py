import sys
import time
import pickle
import multiprocessing
import sqlite3
import numpy as np
import pandas as pd
from desc.gen3_workflow.resource_estimator import SkyMapPolygons, \
    OverlapFinder, extract_coadds, unique_tuples, get_pipetask_resource_funcs, \
    tabulate_pipetask_resources, total_node_hours

# Get skymap convex polygons.
sky_map_file = ('/global/cscratch1/sd/jchiang8/desc/gen3_tests/gen3_repos/'
                'gen3-repo-tiny/skymaps/skyMap/skyMap_DC2_skymaps.pickle')
with open(sky_map_file, 'rb') as fd:
    sky_map = pickle.load(fd)
skymap_polygons = SkyMapPolygons(sky_map, tracts_file='../tracts.pkl')

# Extract first year of WFD minimal dust visits from the opsim db file.
mjd0 = 60218    # 2023-10-01 00:00:00
proposalId = 1
opsim_db_file = '../baseline_v2.0_10yrs.db'
query = f'''select observationId from observations where
            observationStartMJD < {mjd0 + 365}
            and proposalId={proposalId} order by observationId asc'''
with sqlite3.connect(opsim_db_file) as con:
    df = pd.read_sql(query, con)
    visits = list(df['observationId'])

# Simulate the observing cadence and compute overlaps.
overlap_finder = OverlapFinder(opsim_db_file, skymap_polygons,
                               visit_range=(min(visits), max(visits)))

num_visits = len(visits)
batch_size = 500
num_batches = num_visits//batch_size + 1
indexes = [int(_) for _ in np.linspace(96101, num_visits+1, num_batches)]

processes = 64

def write_overlaps(visits, prefix):
    global overlap_finder
    df = overlap_finder.get_overlaps(visits)
    df.to_pickle(f'{prefix}.pickle')

with multiprocessing.Pool(processes=processes) as pool:
    workers = []
    for imin, imax in zip(indexes[:-1], indexes[1:]):
        print(imin, imax)
        sys.stdout.flush()
        prefix = f'overlaps_{imin:05d}_{imax:05d}'
        args = (visits[imin:imax], prefix)
        workers.append(pool.apply_async(write_overlaps, args))
    pool.close()
    pool.join()
    _ = [worker.get() for worker in workers]
