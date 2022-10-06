import os
import glob
from collections import defaultdict
import multiprocessing
import numpy as np
import pandas as pd
import lsst.daf.butler as daf_butler


__all__ = ['get_nImage_stats', 'get_merged_det_stats', 'get_resource_usage',
           'add_nImage_columns', 'add_merged_det_column']


def get_nImage_stats(repo, collection):
    butler = daf_butler.Butler(repo, collections=[collection])
    dstypes = [_ + 'Coadd_nImage' for _ in 'deep goodSeeing'.split()]
    data = defaultdict(list)
    for dstype in dstypes:
        dsrefs = set(butler.registry.queryDatasets(dstype))
        print(dstype, len(dsrefs))
        coadd_type = dstype[:-len('Coadd_nImage')]
        for i, dsref in enumerate(dsrefs):
            print(i)
            image = butler.getDirect(dsref)
            data['coadd_type'].append(coadd_type)
            data['band'].append(dsref.dataId['band'])
            data['tract'].append(dsref.dataId['tract'])
            data['patch'].append(dsref.dataId['patch'])
            data['n_median'].append(np.median(image.array))
            data['n_max'].append(np.max(image.array))
    return pd.DataFrame(data)


def get_merged_det_stats(repo, collection):
    butler = daf_butler.Butler(repo, collections=[collection])
    data = defaultdict(list)
    dstype = 'deepCoadd_mergeDet'
    dsrefs = set(butler.registry.queryDatasets(dstype))
    for dsref in dsrefs:
        cat = butler.getDirect(dsref)
        data['tract'].append(dsref.dataId['tract'])
        data['patch'].append(dsref.dataId['patch'])
        data['n_det'].append(len(cat))
    return pd.DataFrame(data)


def extract_resource_usage(md):
    max_rss_values = []
    for subtask in md.keys():
        for key in md[subtask].keys():
            if 'MaxResidentSetSize' in key:
                # Convert RSS to GB.
                max_rss_values.append(md[subtask][key]/1024**3)
    maxRSS = max(max_rss_values)
    quantum = md['quantum']
    # Convert times to minutes.
    try:
        wall_time = (quantum['endCpuTime'] - quantum['startCpuTime'])/60.
    except KeyError:
        wall_time = None
    try:
        cpu_time = (quantum['endUserTime'] - quantum['startUserTime'])/60.
    except KeyError:
        cpu_time = None
    return maxRSS, wall_time, cpu_time


def fill_data_frames(task, butler, dsrefs):
    data_coadd = defaultdict(list)
    data_visit = defaultdict(list)
    for dsref in dsrefs:
        md = butler.getDirect(dsref)
        try:
            maxRSS, wall_time, cpu_time = extract_resource_usage(md)
        except ValueError:
            pass
        else:
            if 'detector' in dsref.dataId:
                data = data_visit
                data['detector'].append(dsref.dataId['detector'])
                if 'visit' in dsref.dataId:
                    data['visit'].append(dsref.dataId['visit'])
                else:
                    # isr task uses 'exposure' instead of 'visit' (for
                    # DC2 at least), even though it's still the visit.
                    data['visit'].append(dsref.dataId['exposure'])
            elif 'visit' in dsref.dataId:
                # Some focal plane level aggregation tasks pertain to
                # visits, but not individual CCDs.
                data = data_visit
                data['detector'].append(None)
                data['visit'].append(dsref.dataId['visit'])
            elif 'tract' in dsref.dataId:
                data = data_coadd
                data['tract'].append(dsref.dataId['tract'])
                # There are tasks that produce only tract-level outputs.
                data['patch'].append(dsref.dataId.get('patch', None))
            else:
                continue
            data['task'].append(task)
            data['maxRSS (GB)'].append(maxRSS)
            data['wall_time'].append(wall_time)
            data['cpu_time (m)'].append(cpu_time)
            data['band'].append(dsref.dataId.get('band', None))

    df_coadd = pd.DataFrame(data_coadd)
    df_visit = pd.DataFrame(data_visit)
    return df_coadd, df_visit


def get_resource_usage(repo, collections, processes=10, output_name=None,
                       target_dsrefs_size=1000, nmax=None):
    butler = daf_butler.Butler(repo, collections=collections)

    # Find metadata dataset types.
    dstypes = set()
    for collection in collections:
        run_dirs = sorted(glob.glob(os.path.join(repo, collection, '*')))
        for run_dir in run_dirs:
            dstypes.update([os.path.basename(_) for _ in
                            glob.glob(os.path.join(run_dir, '*_metadata'))])

    # Loop over dataset types and extract memory and timing info.
    coadd_dfs = []
    visit_dfs = []
    for i, dstype in enumerate(dstypes):
        task = dstype.split('_')[0]
        dsrefs = list(set(butler.registry.queryDatasets(dstype)))
        if nmax is not None:
            dsrefs = dsrefs[:nmax]
        n_refs = len(dsrefs)
        print(i, task, n_refs, flush=True)
        if n_refs > target_dsrefs_size:
            index = np.linspace(0, n_refs, processes + 1, dtype=int)
            with multiprocessing.Pool(processes=processes) as pool:
                workers = []
                for imin, imax in zip(index[:-1], index[1:]):
                    print('  ', imin, imax, flush=True)
                    workers.append(pool.apply_async(fill_data_frames,
                                                    (task, butler,
                                                     dsrefs[imin:imax])))
                pool.close()
                pool.join()
                for worker in workers:
                    df_coadd, df_visit = worker.get()
                    coadd_dfs.append(df_coadd)
                    visit_dfs.append(df_visit)
        else:
            df_coadd, df_visit = fill_data_frames(task, butler, dsrefs)
            coadd_dfs.append(df_coadd)
            visit_dfs.append(df_visit)

    df_coadd = pd.concat(coadd_dfs)
    df_visit = pd.concat(visit_dfs)

    if output_name is not None:
        df_coadd.to_parquet(f'coadd_resource_usage_{output_name}.parq')
        df_visit.to_parquet(f'visit_resource_usage_{output_name}.parq')

    return df_coadd, df_visit


def add_nImage_columns(df_coadd, df_nImage):
    n_stats = {'deep' : {}, 'goodSeeing' : {}}
    n_ugrizy = {'n_median': defaultdict(lambda : 0),
                'n_max': defaultdict(lambda : 0)}
    for _, row in df_nImage.iterrows():
        n_stats[row.coadd_type][(row.band, row.tract, row.patch)] \
            = (row.n_median, row.n_max)
        if row.coadd_type == 'deep':
            n_ugrizy['n_median'][(row.tract, row.patch)] += row.n_median
            n_ugrizy['n_max'][(row.tract, row.patch)] += row.n_max

    n_median, n_max = [], []
    for _, row in df_coadd.iterrows():
        key = row.band, row.tract, row.patch
        if row.patch is None or np.isnan(row.patch):
            n_median.append(None)
            n_max.append(None)
        elif row.task == 'makeWarp':
            n_median.append(1)
            n_max.append(1)
        elif row.band is None:
            n_median.append(n_ugrizy['n_median'][(row.tract, row.patch)])
            n_max.append(n_ugrizy['n_max'][(row.tract, row.patch)])
        elif row.task == 'templateGen':
            n_median.append(n_stats['goodSeeing'][key][0])
            n_max.append(n_stats['goodSeeing'][key][1])
        else:
            n_median.append(n_stats['deep'][key][0])
            n_max.append(n_stats['deep'][key][1])

    df_coadd['n_median'] = n_median
    df_coadd['n_max'] = n_max

    return df_coadd


def add_merged_det_column(df_coadd, df_merged_det):
    n_det_dict = {(row.tract, row.patch) : row.n_det
                  for _, row in df_merged_det.iterrows()}

    n_det = []
    for _, row in df_coadd.iterrows():
        if row.task == 'deblend':
            key = row.tract, row.patch
            n_det.append(n_det_dict[key])
        else:
            n_det.append(None)
    df_coadd['merged detections'] = n_det
    return df_coadd


if __name__ == "__main__":
    repo = '/global/cfs/cdirs/lsst/production/gen3/DC2/Run2.2i/repo'
    sfp_collection = 'u/descdm/sfp_Y1_3639_visits_part_00_w_2022_38'
    step3_collection = 'u/descdm/coadds_Y1_3639_w_2022_38'
    step4_collection = 'u/descdm/coadds_Y1_3639_w_2022_38_step4'
    step5_collection = 'u/descdm/coadds_Y1_3639_w_2022_38_step5'
    collections = [sfp_collection,
                   step3_collection,
                   step4_collection,
                   step5_collection]

    # nImage info from coadds.
    print("Getting nImage info", flush=True)
    df_nImage = pd.read_parquet('nImage.parq')
    #df_nImage = get_nImage_stats(repo, step3_collection)

    # Merged detection numbers.
    print("Getting merged_det info", flush=True)
    df_merged_det = pd.read_parquet('merged_det.parq')
    #df_merged_det = get_merged_det_stats(repo, step3_collection)

    # Memory and cputime visit- and coadd-level data.
    print("Getting resource usage info", flush=True)
#    df_coadd = pd.read_parquet('coadd.parq')
#    df_visit = pd.read_parquet('visit.parq')
    df_coadd, df_visit = get_resource_usage(repo, collections, nmax=200)

    # Add nImage and merged_det columns to coadd data frame.
    df_coadd = add_nImage_columns(df_coadd, df_nImage)

    # Add merged_det info.
    df_coadd = add_merged_det_column(df_coadd, df_merged_det)
