import sqlalchemy as sql
from sqlalchemy.dialects.sqlite import insert as sqlite_insert
import pandas as pd

import kgenealogic.schema as kg

def increment_batch(engine):
    select_batch = (
        sql.select(kg.meta.c.value)
        .where(kg.meta.c.key=='batch')
    )
    update_batch = (
        sql.update(kg.meta)
        .where(kg.meta.c.key=='batch')
    )
    with engine.connect() as conn:
        batch = int(conn.execute(select_batch).scalar_one())
        batch += 1
        conn.execute(update_batch, [dict(value=str(batch))])
        conn.commit()
    return batch

def add_sources(engine, sources):
    stmt = sqlite_insert(kg.source)
    with engine.connect() as conn:
        conn.execute(stmt, [dict(kit=int(x)) for x in sources])
        conn.commit()

def set_segment_lengths(engine):
    seg_pos = (
        sql.select(
            kg.segment.c["id", "chromosome", "start", "end"],
            sql.func.max(kg.genetmap.c.position).filter(kg.genetmap.c.position <= kg.segment.c.start).label("start1"),
            sql.func.min(kg.genetmap.c.position).filter(kg.genetmap.c.position >= kg.segment.c.start).label("start2"),
            sql.func.max(kg.genetmap.c.position).filter(kg.genetmap.c.position <= kg.segment.c.end).label("end1"),
            sql.func.min(kg.genetmap.c.position).filter(kg.genetmap.c.position >= kg.segment.c.end).label("end2"),
        )
        .join_from(kg.segment, kg.genetmap, kg.segment.c.chromosome==kg.genetmap.c.chromosome)
        .where(kg.segment.c.length.is_(None))
        .group_by(kg.segment.c["id", "chromosome"])
        .cte()
    )

    cm_s1, cm_s2, cm_e1, cm_e2 = kg.genetmap.alias(), kg.genetmap.alias(), kg.genetmap.alias(), kg.genetmap.alias()

    seg_length = (
        sql.select(
            seg_pos.c.id,
            ((cm_e1.c.cm + sql.func.coalesce(((seg_pos.c.end - cm_e1.c.position)/(cm_e2.c.position-cm_e1.c.position))*(cm_e2.c.cm-cm_e1.c.cm),0)) -
            (cm_s1.c.cm + sql.func.coalesce(((seg_pos.c.start - cm_s1.c.position)/(cm_s2.c.position-cm_s1.c.position))*(cm_s2.c.cm-cm_s1.c.cm),0))).label("length")
        )
        .join_from(seg_pos, cm_s1, sql.and_(
            cm_s1.c.chromosome==seg_pos.c.chromosome,
            cm_s1.c.position==sql.func.coalesce(seg_pos.c.start1, seg_pos.c.start2)
        ))
        .join(cm_s2, sql.and_(
            cm_s2.c.chromosome==seg_pos.c.chromosome,
            cm_s2.c.position==sql.func.coalesce(seg_pos.c.start2, seg_pos.c.start1)
        ))
        .join(cm_e1, sql.and_(
            cm_e1.c.chromosome==seg_pos.c.chromosome,
            cm_e1.c.position==sql.func.coalesce(seg_pos.c.end1, seg_pos.c.end2)
        ))
        .join(cm_e2, sql.and_(
            cm_e2.c.chromosome==seg_pos.c.chromosome,
            cm_e2.c.position==sql.func.coalesce(seg_pos.c.end2, seg_pos.c.end1)
        ))
        .subquery()
    )
        
    update_seg_length = (
        sql.update(kg.segment)
        .where(kg.segment.c.id==seg_length.c.id)
        .values(length=seg_length.c.length)
    )
    with engine.connect() as conn:
        conn.execute(update_seg_length)
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

def as_internal_segment(engine, data):
    segments = (
        data[["chromosome", "start", "end"]]
        .groupby(["chromosome", "start", "end"])
        .head(1)
    )
    select_segment_ids = sql.select(kg.segment.c.id.label("segment"), kg.segment.c["chromosome", "start", "end"])
    with engine.connect() as conn:
        segments.to_sql("segment", conn, if_exists="append", index=False)
        conn.commit()
        segment_ids = pd.read_sql(select_segment_ids, conn)

    set_segment_lengths(engine)

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
    add_sources(engine, sources)

    matches = as_internal_segment(engine, matches)

    batch = increment_batch(engine)
    matches['batch'] = batch
    matches = matches[['kit1', 'kit2', 'segment', 'batch']]
    matches = pd.concat([matches, matches.rename(columns=dict(kit1='kit2', kit2='kit1'))],
                        ignore_index=True)

    updated_kits = (
        sql.select(kg.match.c.kit1.distinct().label("kit"))
        .where(kg.match.c.batch==batch)
        .subquery()
    )
    update_sources = (
        sql.update(kg.source)
        .where(kg.source.c.kit==updated_kits.c.kit)
        .values(match=batch)
    )

    with engine.connect() as conn:
        matches.to_sql("match", conn, if_exists="append", index=False)
        conn.execute(update_sources)
        conn.commit()

def import_triangles(engine, triangles):
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
    add_sources(engine, sources)

    triangles = as_internal_segment(engine, triangles)

    batch = increment_batch(engine)
    triangles['batch'] = batch
    triangles = triangles[['kit1', 'kit2', 'kit3', 'segment', 'batch']]
    triangles = pd.concat([
            triangles,
            triangles.rename(columns=dict(kit1='kit2', kit2='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit2='kit3', kit3='kit2')),
            triangles.rename(columns=dict(kit1='kit2', kit2='kit3', kit3='kit1')),
            triangles.rename(columns=dict(kit1='kit3', kit2='kit1', kit3='kit2')),
        ], ignore_index=True)

    updated_kits = (
        sql.select(kg.triangle.c.kit1.distinct().label("kit"))
        .where(kg.triangle.c.batch==batch)
        .subquery()
    )
    update_sources = (
        sql.update(kg.source)
        .where(kg.source.c.kit==updated_kits.c.kit)
        .values(triangle=batch)
    )
    with engine.connect() as conn:
        triangles.to_sql("triangle", conn, if_exists="append", index=False)
        conn.execute(update_sources)
        conn.commit()

__all__ = [import_matches, import_triangles]
