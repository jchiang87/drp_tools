import os
from lsst.daf.butler import Butler, DimensionUniverse
from lsst.pipe.base.graph import QuantumGraph

def tabulate_data_product_sizes(qgraph_file, repo, collection):
    """
    Tabulate the mean sizes of data products listed in a QuantumGraph
    using files in a given repo and collection.

    Parameters
    ----------
    qgraph_file : str
        QuantumGraph file produced by `pipetask qgraph`.
    repo : str
        Path to data repository.
    collection : str
        Collection in repo to use for finding example data products.

    Returns
    -------
    dict(dict(tuple)) Outer dict keyed by task label, inner dicts keyed
    by dataset type with tuple of (mean file size (GB), std file sizes (GB),
    number of files in examples).
    """
    qgraph = QuantumGraph.loadUri(qgraph_file, DimensionUniverse())

    butler = Butler(repo, collections=[collection])
    registry = butler.registry

    # Traverse the QuantumGraph finding the dataset types associated
    # with each task type.
    dstypes = defaultdict(set)
    
    for node in qgraph:
        task = node.taskDef.label
        for dstype in node.quantum.outputs:
            dstypes[task].add(dstype.name)

    dp_inputs = defaultdict(list)
    # Loop over task types and query for each dataset types and
    # compute mean and stdev file sizes for each dataset type.
    data = defaultdict(dict)
    for task, dstypes in dstypes.items():
        for dstype in dstypes:
            file_sizes = [os.stat(butler.getURI(_).path).st_size/1024**3
                          for _ in registry.queryDatasets(dstype)]
            data[task][dstype] = (np.nanmean(file_sizes), np.nanstd(file_sizes),
                                  len(file_sizes))
    return data

root_dir = '/global/cscratch1/sd/jchiang8/desc/gen3_tests'
collection = 'u/jchiang8/drp_3828_24_tiny_sim/20210822T025234Z'
qgraph_file = os.path.join(root_dir, 'w_2021_33', 'submit',
                           f'{collection}/{collection.replace("/", "_")}.qgraph')
repo = os.path.join(root_dir, 'gen3_repos/gen3-3828-y1')

dp_sizes = tabulate_data_product_sizes(qgraph_file, repo, collection)
