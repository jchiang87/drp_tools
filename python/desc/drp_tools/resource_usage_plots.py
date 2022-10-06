from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


__all__ = ['make_visit_resource_usage_plots',
           'make_coadd_resource_usage_plots']


def fit_upper_envelope(x, y):
    df = pd.DataFrame(data=dict(x=x, y=y))
    xvals = sorted(list(set(x)))
    yvals = []
    for xval in xvals:
        yvals.append(np.max(df.query(f'x=={xval}')['y']))
    pars = np.polyfit(xvals, yvals, 1)
    if pars[0] <= 0:
        # Can't have negative slope, so fit a constant line instead.
        pars = (0, np.polyfit(xvals, yvals, 0)[0])
    return pars


def get_percentile_value(xx, percentile=0.95):
    xvals = np.array([_ for _ in xx if _ == _])
    index = int(percentile*len(xvals))
    return sorted(xvals)[index]


def make_visit_resource_usage_plots(df_visit, alpha=1, output_label=None):
    bands = 'ugrizy'
    resource_params = defaultdict(dict)
    for task in sorted(list(set(df_visit['task']))):
        df = df_visit.query(f"task == '{task}'")
        plt.figure(figsize=(8, 4))
        for i, column in enumerate(('cpu_time (m)', 'maxRSS (GB)'), 1):
            plt.subplot(1, 2, i)
            for band in bands:
                my_df = df.query(f"band == '{band}'")
                plt.hist(my_df[column], bins=30, alpha=alpha, label=band)
            pars = (0, get_percentile_value(df[column]))
            resource_params[task][column] = pars
            plt.axvline(pars[1], linestyle='--')
            plt.xlabel(column)
            plt.legend(fontsize='x-small')
            print(task, column, resource_params[task][column])
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(task)
        if output_label is not None:
            outfile = f'{task}_{output_label}.png'
            plt.savefig(outfile)
    return dict(resource_params)


def make_coadd_resource_usage_plots(df_coadd, output_label=None):
    tasks = sorted([_ for _ in set(df_coadd['task']) if
                    ('consolidate' not in _ and 'isolatedStar' not in _)])
    bands = 'ugrizy'

    resource_params = defaultdict(dict)
    for task in tasks:
        df = df_coadd.query(f"task == '{task}'")
        plt.figure(figsize=(8, 4))
        for i, column in enumerate(('cpu_time (m)', 'maxRSS (GB)'), 1):
            plt.subplot(1, 2, i)
            if task == 'deblend':
                x_col = 'merged detections'
                plt.scatter(df[x_col], df[column], s=2)
                plt.xlabel(x_col)
                plt.ylabel(column)
                pars = (0, get_percentile_value(df[column]))
                plt.axhline(pars[1], linestyle='--')
                resource_params[task][column] = pars
            else:
                if len(set(df['band'])) == 1:
                    plt.scatter(df['n_max'], df[column], s=2, label=band)
                    pars = (0, get_percentile_value(df[column]))
                    plt.axhline(pars[1], linestyle='--')
                    resource_params[task][column] = pars
                else:
                    xvals, yvals = [], []
                    for band in bands:
                        my_df = df.query(f"band == '{band}'")
                        if task in ('healSparsePropertyMaps',):
                            plt.hist(my_df[column], bins=30, label=band)
                            xvals.extend(my_df[column])
                        else:
                            xvals.extend(my_df['n_max'])
                            yvals.extend(my_df[column])
                            plt.scatter(my_df['n_max'], my_df[column], s=2,
                                        label=band)
                    if task in ('healSparsePropertyMaps',):
                        pars = (0, get_percentile_value(xvals))
                        plt.axvline(pars[1], linestyle='--')
                        resource_params[task][column] = pars
                    else:
                        plt.ylim(0, 1.05*max(yvals))
                        pars = fit_upper_envelope(xvals, yvals)
                        resource_params[task][column] = tuple(pars)
                        func = np.poly1d(pars)
                        xx = np.linspace(min(xvals), max(xvals), 100)
                        plt.plot(xx, func(xx), linestyle='--')
                plt.legend(fontsize='x-small')
                if task in ('healSparsePropertyMaps',):
                    plt.xlabel(column)
                    plt.ylabel('entries / bin')
                else:
                    plt.xlabel('max(nImage)')
                    plt.ylabel(column)
                plt.tight_layout(rect=(0, 0, 1, 0.95))
            try:
                print(task, column, resource_params[task][column])
            except Exception as eobj:
                print(task, eobj)
        plt.suptitle(task)
        if output_label is not None:
            outfile = f'{task}_{output_label}.png'
            plt.savefig(outfile)
    return dict(resource_params)

if __name__ == '__main__':
    df_coadd = pd.read_parquet('coadd.parq')
    make_coadd_resource_usage_plots(df_coadd)

    df_visit = pd.read_parquet('visit.parq')
    make_visit_resource_usage_plots(df_visit)
