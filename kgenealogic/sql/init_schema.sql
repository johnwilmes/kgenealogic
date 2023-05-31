BEGIN;

/*
    kgenealogic - metadata table
    key: value
        schema_version: version of db schema
        cache_valid: boolean, partition/partition_segment/negative are up to date
*/
DROP TABLE IF EXISTS kgenealogic;
CREATE TABLE kgenealogic (
    key TEXT PRIMARY KEY NOT NULL,
    value TEXT NOT NULL
);
INSERT INTO kgenealogic VALUES
    ('schema_version', '0.1'),
    ('cache_valid', FALSE);

/*
    kit - kit ids and associated metadata
    translation from string kitids to integer ids, for memory efficiency
*/
DROP TABLE IF EXISTS kit;
CREATE TABLE kit (
    id INTEGER PRIMARY KEY NOT NULL,
    kitid TEXT NOT NULL UNIQUE ON CONFLICT IGNORE,
    name TEXT,
    email TEXT,
    sex TEXT
);

/*
    segment - all unique segments that appear in source files
    translate chromosome/start/end to unique segment id
*/
DROP TABLE IF EXISTS segment;
CREATE TABLE segment (
    id INTEGER PRIMARY KEY NOT NULL,
    chromosome TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    length REAL NOT NULL,
    imputed INTEGER NOT NULL DEFAULT FALSE,
    UNIQUE(chromosome, start, end) ON CONFLICT IGNORE
);
CREATE UNIQUE INDEX segment_idx ON segment (id, chromosome, start, end);

/*
    match - pairwise matches on segments
*/
DROP TABLE IF EXISTS match;
CREATE TABLE match (
    segment INTEGER NOT NULL REFERENCES segment (id),
    kit1 INTEGER NOT NULL REFERENCES kit (id),
    kit2 INTEGER NOT NULL REFERENCES kit (id),
    UNIQUE (segment, kit1, kit2) ON CONFLICT IGNORE
);

/*
    triangle - (positive) triangular matches on segments
*/
DROP TABLE IF EXISTS triangle;
CREATE TABLE triangle (
    segment INTEGER NOT NULL REFERENCES segment (id),
    kit1 INTEGER NOT NULL REFERENCES kit (id),
    kit2 INTEGER NOT NULL REFERENCES kit (id),
    kit3 INTEGER NOT NULL REFERENCES kit (id),
    UNIQUE (segment, kit1, kit2, kit3) ON CONFLICT IGNORE
);

COMMIT;
