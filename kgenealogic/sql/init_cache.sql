BEGIN;

DELETE FROM segment WHERE segment.imputed;
UPDATE kgenealogic SET value=FALSE WHERE key='cache_valid';

DROP TABLE IF EXISTS partition;
CREATE TABLE partition (
    id INTEGER PRIMARY KEY NOT NULL,
    chromosome TEXT NOT NULL,
    start INTEGER NOT NULL,
    end INTEGER NOT NULL,
    length REAL,
    mbp REAL GENERATED ALWAYS AS ((end-start)/1e6) STORED,
    UNIQUE(chromosome, start, end) ON CONFLICT IGNORE
);
DROP INDEX IF EXISTS partition_idx;
CREATE UNIQUE INDEX partition_idx ON segment (id, chromosome, start, end);

DROP TABLE IF EXISTS segment_partition;
CREATE TABLE segment_partition (
    segment INTEGER NOT NULL REFERENCES segment (id),
    partition INTEGER NOT NULL REFERENCES partition (id),
    UNIQUE(segment, partition) ON CONFLICT IGNORE
);

/*
    negative - imputed missing triangles
*/
DROP TABLE IF EXISTS negative;
CREATE TABLE negative (
    source INTEGER NOT NULL REFERENCES kit (id),
    target1 INTEGER NOT NULL REFERENCES kit (id),
    target2 INTEGER NOT NULL REFERENCES kit (id),
    segment1 INTEGER NOT NULL REFERENCES segment (id), -- source/target1 positive match
    segment2 INTEGER NOT NULL REFERENCES segment (id), -- source/target2 positive match
    overlap_segment INTEGER NOT NULL REFERENCES segment (id), -- segment1 intersect segment2
    neg_segment INTEGER NOT NULL REFERENCES segment (id), -- subsegment of overlap where target1/target2 don't match
    UNIQUE (source, target1, target2, neg_segment) ON CONFLICT IGNORE
);

COMMIT;
