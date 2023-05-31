BEGIN;
UPDATE kgenealogic SET value=FALSE WHERE key='cache_valid';
COMMIT;
