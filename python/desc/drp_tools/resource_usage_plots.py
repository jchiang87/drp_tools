import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


__all__ = ['make_visit_resource_usage_plots',
           'make_coadd_resource_usage_plots']


def fit_upper_envelope(x, y, degree=1):
    df = pd.DataFrame(data=dict(x=x, y=y))
    xvals = sorted(list(set(x)))
    yvals = []
    for xval in xvals:
        yvals.append(np.max(df.query(f'x=={xval}')['y']))
    pars = np.polyfit(xvals, yvals, degree)
    if pars[0] <= 0:
        # Can't have negative slope, so fit a constant instead.
        pars = np.polyfit(xvals, yvals, degree-1)
    return pars


def make_visit_resource_usage_plots(df_visit):
    bands = 'ugrizy'
    for task in set(df_visit['task']):
        df = df_visit.query(f"task == '{task}'")
        plt.figure(figsize=(8, 4))
        for i, column in enumerate(('cpu_time (min)', 'maxRSS (GB)'), 1):
            plt.subplot(1, 2, i)
            for band in bands:
                my_df = df.query(f"band == '{band}'")
                plt.hist(my_df[column], bins=30, alpha=0.5, label=band)
            plt.xlabel(column)
            plt.legend(fontsize='x-small')
        plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(task)


def make_coadd_resource_usage_plots(df_coadd):
    tasks = sorted([_ for _ in set(df_coadd['task']) if 'consolidate' not in _])
    bands = 'ugrizy'

    for task in tasks:
        df = df_coadd.query(f"task == '{task}'")
        plt.figure(figsize=(8, 4))
        for i, column in enumerate(('cpu_time (min)', 'maxRSS (GB)'), 1):
            plt.subplot(1, 2, i)
            if task == 'deblend':
                x_col = 'merged detections'
                plt.scatter(df[x_col], df[column], s=2)
                plt.xlabel(x_col)
                plt.ylabel(column)
            else:
                if len(set(df['band'])) == 1:
                    plt.scatter(df['n_max'], df[column], s=2, label=band)
                else:
                    xvals, yvals = [], []
                    for band in bands:
                        my_df = df.query(f"band == '{band}'")
                        if task in ('makeWarp', 'healSparsePropertyMaps'):
                            plt.hist(my_df[column], bins=30, alpha=0.4,
                                     label=band)
                        else:
                            xvals.extend(my_df['n_max'])
                            yvals.extend(my_df[column])
                            plt.scatter(my_df['n_max'], my_df[column], s=2,
                                        label=band)
                    if task not in ('makeWarp', 'healSparsePropertyMaps'):
                        plt.ylim(0, 1.05*max(yvals))
                        pars = fit_upper_envelope(xvals, yvals)
                        print(task, pars)
                        func = np.poly1d(pars)
                        xx = np.linspace(min(xvals), max(xvals), 100)
                        plt.plot(xx, func(xx), linestyle='--')
                plt.legend(fontsize='x-small')
                if task in ('makeWarp', 'healSparsePropertyMaps'):
                    plt.xlabel(column)
                    plt.ylabel('entries / bin')
                else:
                    plt.xlabel('max(nImage)')
                    plt.ylabel(column)
                plt.tight_layout(rect=(0, 0, 1, 0.95))
        plt.suptitle(task)

if __name__ == '__main__':
    df_coadd = pd.read_parquet('coadd.parq')
    make_coadd_resource_usage_plots(df_coadd)

    df_visit = pd.read_parquet('visit.parq')
    make_visit_resource_usage_plots(df_visit)
