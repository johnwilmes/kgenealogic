import sqlalchemy as sql
import pandas as pd

from kgenealogic.schema import *

def increment_batch(engine):
    """Increment batch number in meta table and return new value."""
    select_batch = (
        sql.select(meta.c.value)
        .where(meta.c.key=='batch')
    )
    update_batch = (
        sql.update(meta)
        .where(meta.c.key=='batch')
    )
    with engine.connect() as conn:
        batch = int(conn.execute(select_batch).scalar_one())
        batch += 1
        conn.execute(update_batch, [dict(value=str(batch))])
        conn.commit()
    return batch

def add_sources(engine, sources):
    """Add sources to source table, if not already present."""
    with engine.connect() as conn:
        conn.execute(source.insert(), [dict(kit=int(x)) for x in sources])
        conn.commit()

def set_segment_lengths(engine):
    """Compute cM lengths for segments using dnapainter genetic map."""
    seg_pos = (
        sql.select(
            segment.c["id", "chromosome", "start", "end"],
            sql.func.max(genetmap.c.position).filter(genetmap.c.position <= segment.c.start).label("start1"),
            sql.func.min(genetmap.c.position).filter(genetmap.c.position >= segment.c.start).label("start2"),
            sql.func.max(genetmap.c.position).filter(genetmap.c.position <= segment.c.end).label("end1"),
            sql.func.min(genetmap.c.position).filter(genetmap.c.position >= segment.c.end).label("end2"),
        )
        .join_from(segment, genetmap, segment.c.chromosome==genetmap.c.chromosome)
        .where(segment.c.length.is_(None))
        .group_by(segment.c["id", "chromosome"])
        .cte()
    )

    cm_s1, cm_s2, cm_e1, cm_e2 = genetmap.alias(), genetmap.alias(), genetmap.alias(), genetmap.alias()

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
        sql.update(segment)
        .where(segment.c.id==seg_length.c.id)
        .values(length=seg_length.c.length)
    )
    with engine.connect() as conn:
        conn.execute(update_seg_length)
        conn.commit()

def as_internal_kitid(engine, data, kitid_fields):
    """Replace external kitid with internal kit number, inserting into kit table as needed."""
    if len(data)==0:
        return data

    kit_ids = pd.concat([data[c] for c in kitid_fields]).drop_duplicates()
    insert_kit_id = kit.insert().values(kitid=sql.bindparam('kit'))
    select_kit_id = sql.select(kit.c["id", "kitid"])
    with engine.connect() as conn:
        conn.execute(insert_kit_id, [dict(kit=x) for x in kit_ids])
        conn.commit()
        kit_ids = pd.read_sql(select_kit_id, conn)

    kit_ids = kit_ids.set_index('kitid')['id']
    for c in kitid_fields:
        data = data.assign(**{c: data[c].map(kit_ids)})
    return data

def as_internal_segment(engine, data):
    """Replace chromosome/start/end triplet with internal segment id, inserting into segment table
    as needed."""
    if len(data)==0:
        return data.assign(segment=-1)

    seg_df = (
        data[["chromosome", "start", "end"]]
        .groupby(["chromosome", "start", "end"])
        .head(1)
    )
    select_segment_ids = sql.select(segment.c.id.label("segment"), segment.c["chromosome", "start", "end"])
    with engine.connect() as conn:
        conn.execute(segment.insert(), seg_df.to_dict(orient="records"))
        conn.commit()
        segment_ids = pd.read_sql(select_segment_ids, conn)

    set_segment_lengths(engine)

    data = data.merge(segment_ids, on=["chromosome", "start", "end"])
    return data

def update_kit_data(engine, kit_data):
    """Update name/email/sex fields of kit table"""
    insert_kit_data = (
        kit.update()
        .where(kit.c.id==sql.bindparam("kit"))
        .where(kit.c.sex.is_(None))
        .values(
            name=sql.bindparam("name"),
            email=sql.bindparam("email"),
            sex=sql.bindparam("sex"),
        )
    )
    with engine.connect() as conn:
        conn.execute(insert_kit_data, kit_data.to_dict(orient='records'))
        conn.commit()

def import_matches(engine, data):
    """Import pairwise matches to project database from dataframe.

    Args:
        engine: the sqlalchemy engine for the project database
        data: a pandas dataframe with string fields kit1, kit2, chromosome, start, end, name,
        email, sex. Every kit appearing in the kit1 field is assumed to be a "source" for the
        matches
    """
    data = as_internal_kitid(engine, data, ['kit1', 'kit2'])

    kit_data = (
        data[['kit2', 'name', 'email', 'sex']]
        .groupby(['kit2'])
        .head(1)
        .drop_duplicates()
        .rename(columns=dict(kit2='kit'))
    )
    update_kit_data(engine, kit_data)

    sources = data['kit1'].drop_duplicates()
    add_sources(engine, sources)

    data = as_internal_segment(engine, data)

    batch = increment_batch(engine)
    data['batch'] = batch
    data = data[['kit1', 'kit2', 'segment', 'batch']]
    data = pd.concat([data, data.rename(columns=dict(kit1='kit2', kit2='kit1'))],
                        ignore_index=True)

    updated_kits = (
        sql.select(match.c.kit1.distinct().label("kit"))
        .where(match.c.batch==batch)
        .subquery()
    )
    update_sources = (
        sql.update(source)
        .where(source.c.kit==updated_kits.c.kit)
        .values(match=batch)
    )

    with engine.connect() as conn:
        conn.execute(match.insert(), data.to_dict(orient="records"))
        conn.execute(update_sources)
        conn.commit()

def import_triangles(engine, data):
    """Import triangulations to project database from dataframe.

    Args:
        engine: the sqlalchemy engine for the project database
        data: a pandas dataframe with string fields kit1, kit2, kit3, chromosome, start, end, name,
        email. Every kit appearing in the kit1 field is assumed to be a "source" for the triangles
    """
    tri_df = as_internal_kitid(engine, data, ['kit1', 'kit2', 'kit3'])

    kit_data = (
        pd.concat([
            tri_df[['kit2', 'name2', 'email2']].rename(
                columns=dict(kit2='kit', name2='name', email2='email')
            ),
            tri_df[['kit3', 'name3', 'email3']].rename(
                columns=dict(kit3='kit', name3='name', email3='email')
            ),
        ], ignore_index=True)
        .dropna()
        .drop_duplicates()
    )
    kit_data['sex'] = None
    update_kit_data(engine, kit_data)

    sources = tri_df['kit1'].drop_duplicates()
    add_sources(engine, sources)

    tri_df = as_internal_segment(engine, tri_df)

    batch = increment_batch(engine)
    tri_df['batch'] = batch
    tri_df = tri_df[['kit1', 'kit2', 'kit3', 'segment', 'batch']]
    tri_df = pd.concat([
            tri_df,
            tri_df.rename(columns=dict(kit1='kit2', kit2='kit1')),
            tri_df.rename(columns=dict(kit1='kit3', kit3='kit1')),
            tri_df.rename(columns=dict(kit2='kit3', kit3='kit2')),
            tri_df.rename(columns=dict(kit1='kit2', kit2='kit3', kit3='kit1')),
            tri_df.rename(columns=dict(kit1='kit3', kit2='kit1', kit3='kit2')),
        ], ignore_index=True)

    updated_kits = (
        sql.select(triangle.c.kit1.distinct().label("kit"))
        .where(triangle.c.batch==batch)
        .subquery()
    )
    update_sources = (
        sql.update(source)
        .where(source.c.kit==updated_kits.c.kit)
        .values(triangle=batch)
    )
    with engine.connect() as conn:
        conn.execute(triangle.insert(), tri_df.to_dict(orient="records"))
        conn.execute(update_sources)
        conn.commit()

__all__ = [import_matches, import_triangles]
