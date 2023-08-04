import sqlalchemy as sql
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
import pandas as pd

import kgenealogic.schema as kg

def add_sources(engine, sources, has_segments=None, has_triangles=None):
    traits = {}
    if has_segments is not None:
        traits['has_segments'] = has_segments
    if has_triangles is not None:
        traits['has_triangles'] = has_triangles
    stmt = sqlite_insert(kg.source).values(kit=sql.bindparam('kit'), **traits)
    if traits:
        stmt = stmt.on_conflict_do_update(index_elements=['kit'], set_=traits)
    else:
        stmt = stmt.on_conflict_do_nothing(index_elements=['kit'])

    with engine.connect() as conn:
        conn.execute(stmt, [dict(kit=x) for x in sources])
        conn.commit()

def as_internal_kitid(engine, data, kitid_fields):
    kit_ids = pd.concat([data[c] for c in kitid_fields]).drop_duplicates()
    insert_kit_id = kg.kit.insert().values(kitid=sql.bindparam('kit'))
    select_kit_id = sql.select(kg.kit.c["id", "kitid"])
    with engine.connect() as conn:
        conn.execute(insert_kit_id, [dict(kit=x) for x in kit_ids])
        conn.commit()
        kit_ids = pd.read_sql(select_kit_id, conn)

    kit_ids = kit_ids.set_index('kitid')['id']
    for c in kitid_fields:
        data = data.assign(**{c: data[c].map(kit_ids)})
    return data

def as_internal_segment(engine, data, generated=False):
    segments = (
        data[["chromosome", "start", "end"]]
        .groupby(["chromosome", "start", "end"])
        .head(1)
        .assign(generated=generated)
    )
    select_segment_ids = sql.select(kg.segment.c.id.label("segment"), kg.segment.c["chromosome", "start", "end"])
    with engine.connect() as conn:
        segments.to_sql("segment", conn, if_exists="append", index=False)
        conn.commit()
        segment_ids = pd.read_sql(select_segment_ids, conn)

    data = data.merge(segment_ids, on=["chromosome", "start", "end"])
    return data

def update_kit_data(engine, kit_data):
    insert_kit_data = (
        kg.kit.update()
        .where(kg.kit.c.id==sql.bindparam("kit"))
        .where(kg.kit.c.sex.is_(None))
        .values(
            name=sql.bindparam("name"),
            email=sql.bindparam("email"),
            sex=sql.bindparam("sex"),
        )
    )
    with engine.connect() as conn:
        conn.execute(insert_kit_data, kit_data.to_dict(orient='records'))
        conn.commit()

def import_matches(engine, matches):
    """ matches fields: kit1, kit2, chromosome, start, end, name, email, sex"""
    from .cache import invalidate_cache
    invalidate_cache(engine)

    matches = as_internal_kitid(engine, matches, ['kit1', 'kit2'])

    kit_data = (
        matches[['kit2', 'name', 'email', 'sex']]
        .groupby(['kit2'])
        .head(1)
        .drop_duplicates()
        .rename(columns=dict(kit2='kit'))
    )
    update_kit_data(engine, kit_data)

    sources = matches['kit1'].drop_duplicates()
    add_sources(engine, sources, has_segments=True)

    matches = as_internal_segment(engine, matches)

    matches = matches[['kit1', 'kit2', 'segment']]
    matches = pd.concat([matches, matches.rename(columns=dict(kit1='kit2', kit2='kit1'))],
                        ignore_index=True)
    with engine.connect() as conn:
        matches.to_sql("match", conn, if_exists="append", index=False)
        conn.commit()


def import_triangles(engine, triangles):
    from .cache import invalidate_cache
    invalidate_cache(engine)

    triangles = as_internal_kitid(engine, triangles, ['kit1', 'kit2', 'kit3'])

    kit_data = (
        pd.concat([
            triangles[['kit2', 'name2', 'email2']].rename(
                columns=dict(kit2='kit', name2='name', email2='email')
            ),
            triangles[['kit3', 'name3', 'email3']].rename(
                columns=dict(kit3='kit', name3='name', email3='email')
            ),
        ], ignore_index=True)
        .dropna()
        .drop_duplicates()
    )
    kit_data['sex'] = None
    update_kit_data(engine, kit_data)

    sources = triangles['kit1'].drop_duplicates()
    add_sources(engine, sources, has_triangles=True)

    triangles = as_internal_segment(engine, triangles)

    triangles = triangles[['kit1', 'kit2', 'kit3', 'segment']]
    triangles = pd.concat([
            triangles,
            triangles.rename(columns=dict(kit1='kit2', kit2='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit2='kit3', kit3='kit2')),
            triangles.rename(columns=dict(kit1='kit2', kit2='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit2='kit1', kit3='kit2')),
        ], ignore_index=True)
    with engine.connect() as conn:
        triangles.to_sql("triangle", conn, if_exists="append", index=False)
        conn.commit()


__all__ = [import_matches, import_triangles]
