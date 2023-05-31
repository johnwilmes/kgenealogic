import importlib.resources
import pysqlite3 as sql
import pandas as pd

def get_query(name):
    return importlib.resources.files("kgenealogic").joinpath(f"sql/{name}.sql").read_text()

def initialize(db):
    with db:
        db.executescript(get_query("init_schema"))

def is_cache_valid(db):
    is_valid = db.execute("SELECT value FROM kgenealogic WHERE key='cache_valid'").fetchone()[0]
    return bool(int(is_valid))

def import_matches(db, matches):
    """ matches fields: kit1, kit2, chromosome, start, end, length, name, email, sex"""

    kit_ids = pd.concat([matches.kit1, matches.kit2]).drop_duplicates()
    insert_kit_id = "INSERT INTO kit (kitid) VALUES (?)"
    with db:
        db.executescript(get_query("invalidate_cache"))
        db.executemany(insert_kit_id, kit_ids.to_numpy().reshape((-1,1)))

    select_kit_id = "SELECT id, kitid FROM kit"
    kit_ids = pd.read_sql(select_kit_id, db)
    matches = (
        matches
        .merge(kit_ids.rename(columns=dict(kitid='kit1')), on='kit1')
        .drop(columns='kit1')
        .rename(columns=dict(id='kit1'))
        .merge(kit_ids.rename(columns=dict(kitid='kit2')), on='kit2')
        .drop(columns='kit2')
        .rename(columns=dict(id='kit2'))
    )

    kit_data = (
        matches[['kit2', 'name', 'email', 'sex']]
        .groupby(['kit2'])
        .head(1)
        .drop_duplicates()
    )
    insert_kit_data = """
        UPDATE kit SET name=:name, email=:email, sex=:sex
        WHERE id=:kit2 AND sex IS NULL
    """
    with db:
        db.executemany(insert_kit_data, kit_data.to_dict(orient='records'))

    segments = (
        matches[["chromosome", "start", "end", "length"]]
        .groupby(["chromosome", "start", "end"])
        .head(1)
    )
    insert_segments = """
        INSERT INTO segment (chromosome, start, end, length)
        VALUES (:chromosome, :start, :end, :length)
    """
    with db:
        db.executemany(insert_segments, segments.to_dict(orient='records'))
    select_segment_ids = "SELECT id AS segment, chromosome, start, end FROM segment"
    segment_ids = pd.read_sql(select_segment_ids, db)
    matches = matches.merge(segment_ids, on=["chromosome", "start", "end"])

    matches = matches[['kit1', 'kit2', 'segment']]
    matches = pd.concat([matches, matches.rename(columns=dict(kit1='kit2', kit2='kit1'))],
                        ignore_index=True)
    insert_matches = "INSERT INTO match (segment, kit1, kit2) VALUES (:segment, :kit1, :kit2)"
    with db:
        db.executemany(insert_matches, matches.to_dict(orient='records'))

def import_triangles(db, triangles):
    kit_ids = pd.concat([triangles.kit1, triangles.kit2, triangles.kit3]).drop_duplicates()
    insert_kit_id = "INSERT INTO kit (kitid) VALUES (?)"
    with db:
        db.executescript(get_query("invalidate_cache"))
        db.executemany(insert_kit_id, kit_ids.to_numpy().reshape((-1,1)))

    select_kit_id = "SELECT id, kitid FROM kit"
    kit_ids = pd.read_sql(select_kit_id, db)
    triangles = (
        triangles
        .merge(kit_ids.rename(columns=dict(kitid='kit1')), on='kit1')
        .drop(columns='kit1')
        .rename(columns=dict(id='kit1'))
        .merge(kit_ids.rename(columns=dict(kitid='kit2')), on='kit2')
        .drop(columns='kit2')
        .rename(columns=dict(id='kit2'))
        .merge(kit_ids.rename(columns=dict(kitid='kit3')), on='kit3')
        .drop(columns='kit3')
        .rename(columns=dict(id='kit3'))
    )

    kit_data = (
        pd.concat([
            triangles[['kit1', 'name1', 'email1']].rename(
                columns=dict(kit1='kit', name1='name', email1='email')
            ),
            triangles[['kit2', 'name2', 'email2']].rename(
                columns=dict(kit2='kit', name2='name', email2='email')
            ),
        ], ignore_index=True)
        .dropna()
    )
    insert_kit_data = """
        UPDATE kit SET name=:name, email=:email
        WHERE id=:kit AND name IS NULL
    """
    with db:
        db.executemany(insert_kit_data, kit_data.to_dict(orient='records'))

    segments = (
        triangles[["chromosome", "start", "end", "length"]]
        .groupby(["chromosome", "start", "end"])
        .head(1)
    )
    insert_segments = """
        INSERT INTO segment (chromosome, start, end, length)
        VALUES (:chromosome, :start, :end, :length)
    """
    with db:
        db.executemany(insert_segments, segments.to_dict(orient='records'))
    select_segment_ids = "SELECT id AS segment, chromosome, start, end FROM segment"
    segment_ids = pd.read_sql(select_segment_ids, db)
    triangles = triangles.merge(segment_ids, on=["chromosome", "start", "end"])

    triangles = triangles[['kit1', 'kit2', 'kit3', 'segment']]
    triangles = pd.concat([
            triangles,
            triangles.rename(columns=dict(kit1='kit2', kit2='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit2='kit3', kit3='kit2')),
            triangles.rename(columns=dict(kit1='kit2', kit2='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit2='kit1', kit3='kit2')),
        ], ignore_index=True)
    insert_triangles ="""
        INSERT INTO triangle (segment, kit1, kit2, kit3)
        VALUES (:segment, :kit1, :kit2, :kit3)
    """
    with db:
        db.executemany(insert_triangles, triangles.to_dict(orient='records'))
