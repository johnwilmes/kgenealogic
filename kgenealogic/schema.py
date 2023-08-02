import sqlalchemy as sql
import pandas as pd

SCHEMA_VERSION = 0.1
metadata = sql.MetaData()

kgenealogic = sql.Table(
    "kgenealogic",
    metadata,
    sql.Column("key", sql.String, nullable=False, primary_key=True),
    sql.Column("value", sql.String, nullable=False),
)

kit = sql.Table(
    "kit",
    metadata,
    sql.Column("id", sql.Integer, nullable=False, primary_key=True),
    sql.Column("kitid", sql.String, nullable=False, unique=True, sqlite_on_conflict_unique='IGNORE'),
    sql.Column("name", sql.String),
    sql.Column("email", sql.String),
    sql.Column("sex", sql.String),
)

source = sql.Table(
    "source",
    metadata,
    sql.Column("kit", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, primary_key=True),
    sql.Column("has_segments", sql.Boolean, nullable=False, default=False),
    sql.Column("has_triangles", sql.Boolean, nullable=False, default=False),
    sql.Column("has_negative", sql.Boolean, nullable=False, default=False),
)

segment = sql.Table(
    "segment",
    metadata,
    sql.Column("id", sql.Integer, nullable=False, primary_key=True),
    sql.Column("chromosome", sql.String, nullable=False, index=True),
    sql.Column("start", sql.Integer, nullable=False),
    sql.Column("end", sql.Integer, nullable=False),
    sql.Column("length", sql.Float, index=True),
    sql.Column("generated", sql.Boolean, nullable=False),
    sql.UniqueConstraint("chromosome", "start", "end", sqlite_on_conflict='IGNORE'),
)

match = sql.Table(
    "match",
    metadata,
    sql.Column("segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.Column("kit1", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit2", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.UniqueConstraint("segment", "kit1", "kit2", sqlite_on_conflict='IGNORE'),
)

triangle = sql.Table(
    "triangle",
    metadata,
    sql.Column("segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.Column("kit1", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit2", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit3", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.UniqueConstraint("segment", "kit1", "kit2", "kit3", sqlite_on_conflict='IGNORE'),
)

overlap = sql.Table(
    "overlap",
    metadata,
    sql.Column("id", sql.Integer, nullable=False, primary_key=True),
    sql.Column("source", sql.Integer, sql.ForeignKey("source.kit"), nullable=False, index=True),
    sql.Column("target1", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("target2", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.UniqueConstraint("source", "target1", "target2", "segment", sqlite_on_conflict='IGNORE'),
)

negative = sql.Table(
    "negative",
    metadata,
    sql.Column("overlap", sql.Integer, sql.ForeignKey("overlap.id"), nullable=False, index=True),
    sql.Column("neg_segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.UniqueConstraint("overlap", "neg_segment", sqlite_on_conflict='IGNORE'),
)

def initialize(engine):
    metadata.drop_all(engine)
    metadata.create_all(engine)
    kg_meta = pd.DataFrame([
        dict(key='schema_version', value=str(SCHEMA_VERSION)),
        dict(key='cache_valid', value=str('N')),
    ])
    with engine.connect() as conn:
        kg_meta.to_sql("kgenealogic", conn, if_exists="append", index=False)
        conn.commit()
