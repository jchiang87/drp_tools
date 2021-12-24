"""
Functions to fill visit tract overlap tables for organizing an
image processing campaign.
"""
import os
import pickle
import sqlite3
import numpy as np
import pandas as pd
import lsst.geom

__all__ = ['fill_visit_table', 'fill_tract_table', 'find_visit_tract_overlaps',
           'fill_overlap_table', 'update_visit_table']

# DC2 tracts, 151 total
DC2_TRACTS = []
DC2_TRACTS.extend(range(2897, 2909))
DC2_TRACTS.extend(range(3074, 3087))
DC2_TRACTS.extend(range(3256, 3269))
DC2_TRACTS.extend(range(3442, 3454))
DC2_TRACTS.extend(range(3631, 3644))
DC2_TRACTS.extend(range(3825, 3838))
DC2_TRACTS.extend(range(4023, 4035))
DC2_TRACTS.extend(range(4224, 4237))
DC2_TRACTS.extend(range(4429, 4441))
DC2_TRACTS.extend(range(4636, 4649))
DC2_TRACTS.extend(range(4848, 4861))
DC2_TRACTS.extend(range(5062, 5075))

VISIT_TABLE = 'Visit'
TRACT_TABLE = 'Tract'
OVERLAP_TABLE = 'Overlap'


def max_tract_radius(tracts=DC2_TRACTS, skymap_file='/home/DC2/skyMap.pickle'):
    """
    Maximum distance (in degrees) from vertex to tract center for the
    list of tracts supplied.  For DC2 tracts, this is 1.1 degrees.
    """
    with open(skymap_file, 'rb') as fobj:
        skymap = pickle.load(fobj)
    max_radius = None
    for tract_id in DC2_TRACTS:
        tract = skymap[tract_id]
        center = tract.ctr_coord
        for vertex in tract.vertex_list:
            sep = center.separation(vertex)
            if max_radius is None or sep > max_radius:
                max_radius = sep
    return max_radius.asDegrees()


def fill_tract_table(db_file, skymap_file='/home/DC2/skyMap.pickle',
                     tract_table=TRACT_TABLE, tract_list=DC2_TRACTS):
    """
    Fill the tract table with the tract ids and tract centers for the
    provided tract_list.
    """
    with open(skymap_file, 'rb') as fobj:
        skymap = pickle.load(fobj)

    create_tract_table = (f'create table {tract_table} '
                          '(id INTEGER, ra REAL, dec REAL)')
    with sqlite3.connect(db_file) as con:
        cursor = con.cursor()
        cursor.execute(create_tract_table)
        con.commit()
        values = []
        for tract_id in tract_list:
            tract = skymap[tract_id]
            coord = tract.getCtrCoord()
            values.append((tract_id, coord.getLongitude().asDegrees(),
                           coord.getLatitude().asDegrees()))
        cursor.executemany(f'insert into {tract_table} values (?, ?, ?)',
                           values)
        con.commit()


def fill_visit_table(db_file, opsim_db='/home/DC2/minion_1016_desc_dithered_v4_trimmed.db', visit_table=VISIT_TABLE):
    """
    Fill the visit table with entries from the opsim db and
    provide a column for the nearest tract to use for partitioning
    the visits for processing.
    """
    visit_table_sql = (f'create table {visit_table} '
                       '(id INTEGER, ra REAL, dec REAL, band TEXT, '
                       'nearest_tract INTEGER, survey_id INT, mjd REAL)')
    with sqlite3.connect(opsim_db) as con:
        df = pd.read_sql('''select obsHistID, descDitheredRA, descDitheredDec,
                         propID, expMJD, filter from Summary''', con)
    with sqlite3.connect(db_file) as con:
        cursor = con.cursor()
        cursor.execute(visit_table_sql)
        con.commit()
        values = []
        for i, row in df.iterrows():
            values.append((row['obsHistID'], row['descDitheredRA']*180/np.pi,
                           row['descDitheredDec']*180/np.pi, row['filter'],
                           0, row['propID'], row['expMJD']))
        cursor.executemany(f'insert into {visit_table} values '
                           '(?, ?, ?, ?, ?, ?, ?)', values)
        con.commit()


def find_visit_tract_overlaps(ra0, dec0, df_tracts, max_sep=3.15):
    """
    Using the maximum separation between an overlapping visit and
    tract, find all tracts overlapping the visit whose center
    coordinate is provided. Also find the closest tract to that center
    coordinate among all tracts in df_tracts.
    """
    dec_min, dec_max = dec0 - max_sep, dec0 + max_sep
    ra_delta = max_sep/np.cos(dec0*np.pi/180.)
    ra_min, ra_max = ra0 - ra_delta, ra0 + ra_delta
    df = df_tracts.query(f'{ra_min} <= ra and ra <= {ra_max} and '
                         f'{dec_min} <= dec and dec <= {dec_max}')
    def make_sphere_point(ra, dec):
        ra_angle = lsst.geom.Angle(ra, lsst.geom.degrees)
        dec_angle = lsst.geom.Angle(dec, lsst.geom.degrees)
        return lsst.geom.SpherePoint(ra_angle, dec_angle)
    visit_center = make_sphere_point(ra0, dec0)
    min_sep = None
    tract_overlaps = set()
    for i, row in df.iterrows():
        tract_center = make_sphere_point(row.ra, row.dec)
        sep = visit_center.separation(tract_center).asDegrees()
        if min_sep is None or sep < min_sep:
            min_sep = sep
            closest_tract = row.id
        if sep <= max_sep:
            tract_overlaps.add(row.id)
    if len(df) == 0:
        for i, row in df_tracts.iterrows():
            tract_center = make_sphere_point(row.ra, row.dec)
            sep = visit_center.separation(tract_center).asDegrees()
            if min_sep is None or sep < min_sep:
                min_sep = sep
                closest_tract = row.id
    return tract_overlaps, closest_tract


def fill_overlap_table(db_file, overlap_table=OVERLAP_TABLE,
                       visit_table=VISIT_TABLE, tract_table=TRACT_TABLE,
                       max_sep=3.15):
    """
    Fill the Overlap table which lists all of the potential overlapping
    visit-tract pairs.  Also provide the dict of closest tracts for each
    visit.
    """
    overlap_table_sql = (f'create table if not exists {overlap_table} '
                         '(id INTEGER, tract INTEGER, visit INTEGER)')
    with sqlite3.connect(db_file) as con:
        df_visits = pd.read_sql(f'select * from {visit_table}', con)
        df_tracts = pd.read_sql(f'select * from {tract_table}', con)
        cursor = con.cursor()
        cursor.execute(overlap_table_sql)
        con.commit()
        values = []
        id_ = 0
        closest_tracts = dict()
        num_visits = len(df_visits)
        for i, visit in df_visits.iterrows():
            if i % (num_visits//5) == 0:
                print('!', end='', flush=True)
            elif i % (num_visits//20) == 0:
                print('.', end='', flush=True)
            tract_overlaps, closest_tract \
                = find_visit_tract_overlaps(visit['ra'], visit['dec'],
                                            df_tracts, max_sep=max_sep)
            closest_tracts[visit['id']] = closest_tract
            for tract in tract_overlaps:
                values.append((id_, tract, visit['id']))
                id_ += 1
        cursor.executemany((f'insert into {overlap_table} values '
                           '(?, ?, ?)'), values)
    return closest_tracts


def update_visit_table(db_file, closest_tracts, visit_table=VISIT_TABLE):
    """
    Update the visit table with the closest tracts.
    """
    with sqlite3.connect(db_file) as con:
        sql = f'update {visit_table} set nearest_tract=? where id=?'
        values = [(tract, visit) for visit, tract in closest_tracts.items()]
        con.cursor().executemany(sql, values)
        con.commit()


if __name__ == '__main__':
    db_file = 'drp_tables.db'
    fill_visit_table(db_file)
    fill_tract_table(db_file)
    closest_tracts = fill_overlap_table(db_file)
    with open('closest_tracts.pickle', 'wb') as fobj:
        pickle.dump(closest_tracts, fobj)
    with open('closest_tracts.pickle', 'rb') as fobj:
        closest_tracts = pickle.load(fobj)
    update_visit_table(db_file, closest_tracts)
