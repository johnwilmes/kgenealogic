import sqlalchemy as sql
import pandas as pd
import importlib.resources as resources

SCHEMA_VERSION = "0.2"
GENETMAP_PATH='genetmap.csv'
GENETMAP_DTYPE={'chromosome': str, 'position': int, 'cm': float}

@sql.event.listens_for(sql.engine.Engine, "connect")
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()

metadata = sql.MetaData()

meta = sql.Table(
    "meta",
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

genetmap = sql.Table(
    "genetmap",
    metadata,
    sql.Column("chromosome", sql.String, nullable=False, index=True),
    sql.Column("position", sql.Integer, nullable=False, index=True),
    sql.Column("cm", sql.Float, nullable=False),
    sql.UniqueConstraint("chromosome", "position"),
)

source = sql.Table(
    "source",
    metadata,
    sql.Column("kit", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, primary_key=True,
               sqlite_on_conflict_primary_key="IGNORE"),
    sql.Column("match", sql.Integer), # batch number of most recent matches
    sql.Column("triangle", sql.Integer), # batch number of most recent triangles
    sql.Column("negative", sql.Integer), # negative triangles incorporate data from this batch
)

segment = sql.Table(
    "segment",
    metadata,
    sql.Column("id", sql.Integer, nullable=False, primary_key=True),
    sql.Column("chromosome", sql.String, nullable=False, index=True), # chromosome 23 is "X"
    sql.Column("start", sql.Integer, nullable=False),
    sql.Column("end", sql.Integer, nullable=False),
    sql.Column("length", sql.Float, index=True),
    sql.UniqueConstraint("chromosome", "start", "end", sqlite_on_conflict='IGNORE'),
)

# pairwise matches are included with both permutations of (kit1, kit2)
match = sql.Table(
    "match",
    metadata,
    sql.Column("segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.Column("kit1", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit2", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("batch", sql.Integer, nullable=False, index=True), # when it was added
    sql.UniqueConstraint("segment", "kit1", "kit2", sqlite_on_conflict='IGNORE'),
)

# triangles are included with all 6 permutations of (kit1, kit2, kit3)
triangle = sql.Table(
    "triangle",
    metadata,
    sql.Column("segment", sql.Integer, sql.ForeignKey("segment.id"), nullable=False, index=True),
    sql.Column("kit1", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit2", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("kit3", sql.Integer, sql.ForeignKey("kit.id"), nullable=False, index=True),
    sql.Column("batch", sql.Integer, nullable=False, index=True), # when it was added
    sql.UniqueConstraint("segment", "kit1", "kit2", "kit3", sqlite_on_conflict='IGNORE'),
)

# gives segments where (source, target1) and (source, target2) have overlapping matches
# used for finding negative triangles
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
    sql.Column("overlap", sql.Integer, sql.ForeignKey("overlap.id", ondelete="CASCADE"),
               nullable=False, index=True),
    sql.Column("neg_segment", sql.Integer, sql.ForeignKey("segment.id", ondelete="CASCADE"),
               nullable=False, index=True),
    sql.UniqueConstraint("overlap", "neg_segment", sqlite_on_conflict='IGNORE'),
)

def initialize(engine):
    metadata.drop_all(engine)
    metadata.create_all(engine)

    kg_meta_init = [
        dict(key='schema_version', value=str(SCHEMA_VERSION)),
        dict(key='batch', value=str(0)),
    ]

    import kgenealogic as kg
    genetmap_df = pd.read_csv(resources.files(kg).joinpath(GENETMAP_PATH), dtype=GENETMAP_DTYPE)

    with engine.connect() as conn:
        conn.execute(meta.insert(), kg_meta_init)
        conn.execute(genetmap.insert(), genetmap_df.to_dict(orient="records"))
        conn.commit()
