-- Jump RV queries for the DR4 study, written against Jump's schema page (2026-11-04) and
-- first run on the pilot systems that day. Paste each query into Jump's SQL Explorer
-- (/sql) and export the result as CSV into rv_data/jump_raw/ (A: jump_matches_<date>.csv,
-- B: jump_rvs_<date>.csv), then run rv_data/ingest_jump.py. To change the targets,
-- replace the VALUES list in both queries with rv_data/targets_values.sql, which
-- rv_data/build_targets.py writes.


-- =====================================================================================
-- A. Name matching: one row per (target, Jump star)
-- =====================================================================================
-- A target with an empty star_id matched nothing; a target with two rows matched two
-- stars. Names are compared in upper case without spaces, hyphens, or leading zeros.
-- Jump names come from core_star, core_staralias, the SIMBAD cross-match, and the TIC
-- and TOI numbers in core_tic and core_toi.

WITH targets(target, aliases) AS (
  VALUES
  ('TOI-1694', 'TOI-1694; TIC 396740648; Gaia DR3 1104269596843209472; 2MASS J06305955+6621384; TYC 4108-1434-1'),
  ('HD 111232', 'HD 111232; HIP 62534; TIC 329649885; Gaia DR3 5855730584310531200; 2MASS J12485177-6825304; TYC 9241-249-1'),
  ('HD 191939', 'HD 191939; HIP 99175; TOI-1339; TIC 269701147; Gaia DR3 2248126315275354496; 2MASS J20080574+6651019; TYC 4244-964-1'),
  ('TOI-1736', 'TOI-1736; TIC 408618999; Gaia DR3 541725187117160960; 2MASS J02534440+6906050; TYC 4313-1054-1'),
  ('HD 43691', 'HD 43691; HIP 30057; TIC 258444495; Gaia DR3 959971546241608576; 2MASS J06193467+4105321; TYC 2930-2105-1')
),
t AS (
  SELECT target, trim(a) AS alias,
         regexp_replace(regexp_replace(upper(trim(a)), '(^|[^0-9])0+([0-9])', '\1\2', 'g'),
                        '[ -]', '', 'g') AS k
  FROM targets, regexp_split_to_table(aliases, ';') AS a
),
names AS (
  SELECT id AS star_id, name AS nm FROM core_star
  UNION SELECT star_id, name FROM core_staralias
  UNION SELECT star_id, simbad_name FROM core_simbadcrossmatch
  UNION SELECT star_id, resolvable_name FROM core_simbadcrossmatch
  UNION SELECT star_id, CASE WHEN tic ~ '^[0-9]+$' THEN 'TIC ' || tic ELSE tic END FROM core_tic
  UNION SELECT star_id, 'TOI-' || substring(toi FROM '[0-9]+') FROM core_tic WHERE toi ~ '[0-9]'
  UNION SELECT star_id, 'TOI-' || substring(toi_cand FROM '[0-9]+') FROM core_toi WHERE toi_cand ~ '[0-9]'
),
n AS (
  -- Names of a merged star point to the star it was merged into
  SELECT COALESCE(s.merged_into_id, s.id) AS star_id, names.nm,
         regexp_replace(regexp_replace(upper(names.nm), '(^|[^0-9])0+([0-9])', '\1\2', 'g'),
                        '[ -]', '', 'g') AS k
  FROM names JOIN core_star s ON s.id = names.star_id
  WHERE names.nm IS NOT NULL
)
SELECT t.target, n.star_id, s.name AS jump_star_name,
       string_agg(DISTINCT t.alias, '; ') AS matched_alias,
       string_agg(DISTINCT n.nm, '; ') AS matched_jump_names,
       s.ra, s.dec, s.vmag,
       s.num_observations_hires, s.num_observations_apf, s.num_observations_kpf,
       s.num_nights_imported, s.has_external_rvs
FROM t
LEFT JOIN n ON n.k = t.k
LEFT JOIN core_star s ON s.id = n.star_id
GROUP BY t.target, n.star_id, s.name, s.ra, s.dec, s.vmag, s.num_observations_hires, s.num_observations_apf,
         s.num_observations_kpf, s.num_nights_imported, s.has_external_rvs
ORDER BY t.target, s.name;


-- =====================================================================================
-- B. RVs: one row per measurement, from every RV table in Jump
-- =====================================================================================
-- Sources, each tagged in `database` (the two-tier provenance's first tier; the
-- `publication` column is the second):
--   HIRES and APF    core_pipelinervmb / core_pipelinervapf, latest reduction of each
--                    exposure (by ingest_time); iodine-free exposures (templates) left out
--   (KPF is left out on purpose, decided 2026-11-04: the instrument isn't yet stable
--   over the long term. Its RVs are in core_pipelinervkpf, which the SQL Explorer
--   doesn't show anyway.)
--   NEID             core_neidrv, and core_proprietaryneidrv (with its status and embargo)
--   literature       core_literaturerv (with source and bibcode)
--   imports          core_externalrv (DACE, HARPS-RVBank, OHP, VizieR, uploads), with the
--                    dataset's source, time scale, corrections, unit, bibcode, and status
-- utctime is the observation table's UTC time, which for HIRES (after 2004) and APF
-- appears to be the exposure start: it precedes the BJD-implied midpoint by up to about
-- half the exposure. jd_kbcvel is the HIRES midpoint (JD, UTC).
-- HIRES and APF time_bjd are BJD_TDB (checked 2026-11-04: HIRES to < 0.01 s against
-- jd_kbcvel; APF to about 2 s against utctime + exptime/2, against 69 s for UTC).
-- Quality flags are returned, not applied. Every dataset row is returned, including
-- retired or superseded ones, with their status in quality_flag, so duplicates across
-- datasets can be resolved afterwards.

WITH targets(target, aliases) AS (
  VALUES
  ('TOI-1694', 'TOI-1694; TIC 396740648; Gaia DR3 1104269596843209472; 2MASS J06305955+6621384; TYC 4108-1434-1'),
  ('HD 111232', 'HD 111232; HIP 62534; TIC 329649885; Gaia DR3 5855730584310531200; 2MASS J12485177-6825304; TYC 9241-249-1'),
  ('HD 191939', 'HD 191939; HIP 99175; TOI-1339; TIC 269701147; Gaia DR3 2248126315275354496; 2MASS J20080574+6651019; TYC 4244-964-1'),
  ('TOI-1736', 'TOI-1736; TIC 408618999; Gaia DR3 541725187117160960; 2MASS J02534440+6906050; TYC 4313-1054-1'),
  ('HD 43691', 'HD 43691; HIP 30057; TIC 258444495; Gaia DR3 959971546241608576; 2MASS J06193467+4105321; TYC 2930-2105-1')
),
t AS (
  SELECT target, trim(a) AS alias,
         regexp_replace(regexp_replace(upper(trim(a)), '(^|[^0-9])0+([0-9])', '\1\2', 'g'),
                        '[ -]', '', 'g') AS k
  FROM targets, regexp_split_to_table(aliases, ';') AS a
),
names AS (
  SELECT id AS star_id, name AS nm FROM core_star
  UNION SELECT star_id, name FROM core_staralias
  UNION SELECT star_id, simbad_name FROM core_simbadcrossmatch
  UNION SELECT star_id, resolvable_name FROM core_simbadcrossmatch
  UNION SELECT star_id, CASE WHEN tic ~ '^[0-9]+$' THEN 'TIC ' || tic ELSE tic END FROM core_tic
  UNION SELECT star_id, 'TOI-' || substring(toi FROM '[0-9]+') FROM core_tic WHERE toi ~ '[0-9]'
  UNION SELECT star_id, 'TOI-' || substring(toi_cand FROM '[0-9]+') FROM core_toi WHERE toi_cand ~ '[0-9]'
),
n AS (
  SELECT COALESCE(s.merged_into_id, s.id) AS star_id, names.nm,
         regexp_replace(regexp_replace(upper(names.nm), '(^|[^0-9])0+([0-9])', '\1\2', 'g'),
                        '[ -]', '', 'g') AS k
  FROM names JOIN core_star s ON s.id = names.star_id
  WHERE names.nm IS NOT NULL
),
m AS (
  SELECT t.target, n.star_id, string_agg(DISTINCT t.alias, '; ') AS matched_alias
  FROM t JOIN n ON n.k = t.k
  GROUP BY t.target, n.star_id
),
-- Every row of a star-linked table, mapped through merges to the matched star
star_map AS (
  SELECT m.target, m.matched_alias, m.star_id, s.id AS any_star_id
  FROM m JOIN core_star s ON COALESCE(s.merged_into_id, s.id) = m.star_id
),
-- Observations of the matched stars (match columns renamed so they can't collide)
obs AS (
  SELECT sm.target AS m_target, sm.matched_alias AS m_alias, sm.star_id AS m_star_id, o.*
  FROM core_observation o JOIN star_map sm ON sm.any_star_id = o.star_id
),
legs AS (
  -- HIRES (instrument hires_k before the August 2004 upgrade, hires_j after)
  SELECT o.m_target AS target, o.m_alias AS matched_alias, o.m_star_id AS star_id,
         'Jump: CPS pipeline'::text AS database, o.observation_id::text AS obs_id,
         o.instrument::text AS instrument, k.bjd AS time_bjd, 'unknown'::text AS time_scale,
         r.mnvel AS rv_ms, r.errvel AS rv_err_ms,
         concat_ws('; ', 'vst_file=' || r.vst_file, 'ingest_time=' || r.ingest_time) AS reduction,
         NULL::text AS corrections,
         concat_ws('; ', 'isjunk=' || o.isjunk, 'mdchi=' || r.mdchi) AS quality_flag,
         concat_ws('; ', 'cts=' || r.cts, 'hdr_counts=' || o.hdr_counts, 'hdr_snr=' || o.hdr_snr) AS counts_or_snr,
         sv.svalue AS s_index, o.program_id AS program_id, o.raw_progid::text AS program_code,
         NULL::text AS publication,
         concat_ws('; ', 'program_locked=' || o.program_locked, 'embargo_locked=' || o.embargo_locked,
                   'embargo_expiration=' || o.embargo_expiration) AS access,
         o.utctime::text AS utctime, k.jd AS jd_kbcvel, 'm/s'::text AS rv_unit,
         o.hdr_exposure_time AS exptime_s
  FROM obs o
  JOIN LATERAL (SELECT * FROM core_pipelinervmb p WHERE p.observation_id = o.id
                ORDER BY p.ingest_time DESC LIMIT 1) r ON true
  LEFT JOIN LATERAL (SELECT bjd, jd FROM core_kbcvel b WHERE b.observation_id = o.id
                     ORDER BY b.id DESC LIMIT 1) k ON true
  LEFT JOIN LATERAL (SELECT svalue FROM core_sval v WHERE v.observation_id = o.id
                     ORDER BY v.id DESC LIMIT 1) sv ON true
  WHERE o.iodine_in IS NOT FALSE

  UNION ALL
  -- APF
  SELECT o.m_target, o.m_alias, o.m_star_id, 'Jump: CPS pipeline', o.observation_id,
         o.instrument, k.bjd, 'unknown', r.mnvel, r.errvel,
         concat_ws('; ', 'vst_file=' || r.vst_file, 'ingest_time=' || r.ingest_time), NULL,
         concat_ws('; ', 'isjunk=' || o.isjunk, 'mdchi=' || r.mdchi),
         concat_ws('; ', 'cts=' || r.cts, 'hdr_counts=' || o.hdr_counts, 'hdr_snr=' || o.hdr_snr),
         sv.svalue, o.program_id, o.raw_progid, NULL,
         concat_ws('; ', 'program_locked=' || o.program_locked, 'embargo_locked=' || o.embargo_locked,
                   'embargo_expiration=' || o.embargo_expiration),
         o.utctime::text, NULL::float8, 'm/s', o.hdr_exposure_time
  FROM obs o
  JOIN LATERAL (SELECT * FROM core_pipelinervapf p WHERE p.observation_id = o.id
                ORDER BY p.ingest_time DESC LIMIT 1) r ON true
  LEFT JOIN LATERAL (SELECT bjd FROM core_abcvel b WHERE b.observation_id = o.id
                     ORDER BY b.id DESC LIMIT 1) k ON true
  LEFT JOIN LATERAL (SELECT svalue FROM core_svalapf v WHERE v.observation_id = o.id
                     ORDER BY v.id DESC LIMIT 1) sv ON true
  WHERE o.iodine_in IS NOT FALSE

  UNION ALL
  -- NEID linked to observations (time stored as `jd`, scale not recorded)
  SELECT o.m_target, o.m_alias, o.m_star_id, 'Jump: NEID', o.observation_id,
         o.instrument, r.jd, 'unknown (jd column)', r.mnvel, r.errvel,
         'swversion=' || r.swversion, NULL,
         'isjunk=' || o.isjunk, 'snr=' || r.snr,
         NULL::float8, o.program_id, o.raw_progid, NULL,
         concat_ws('; ', 'program_locked=' || o.program_locked, 'embargo_locked=' || o.embargo_locked,
                   'embargo_expiration=' || o.embargo_expiration),
         o.utctime::text, NULL::float8, NULL, o.hdr_exposure_time
  FROM obs o
  JOIN core_neidrv r ON r.observation_id = o.id

  UNION ALL
  -- Proprietary NEID archive rows
  SELECT sm.target, sm.matched_alias, sm.star_id, 'Jump: proprietary NEID',
         'neid:' || r.archive_row_id, 'neid', r.epoch_bjd, 'unknown', r.rv, r.rv_err,
         'drp_version=' || r.drp_version, NULL,
         concat_ws('; ', 'status=' || r.status, 'caution=' || r.caution), NULL,
         NULL::float8, r.program_id, r.program_code, NULL,
         'embargo_expiration=' || r.embargo_expiration,
         NULL, NULL::float8, NULL, NULL::float8
  FROM core_proprietaryneidrv r JOIN star_map sm ON sm.any_star_id = r.star_id

  UNION ALL
  -- Literature RVs stored in Jump (time stored as `jd`, scale not recorded)
  SELECT sm.target, sm.matched_alias, sm.star_id,
         'Jump: literature RVs (' || COALESCE(r.source, 'source not recorded') || ')',
         'lit:' || r.id, r.instrument, r.jd, 'unknown (jd column)', r.mnvel, r.errvel,
         'inst_code=' || r.inst_code, NULL, NULL, 'snr=' || r.snr,
         NULL::float8, NULL::bigint, NULL, r.bibcode, NULL,
         NULL, NULL::float8, NULL, NULL::float8
  FROM core_literaturerv r JOIN star_map sm ON sm.any_star_id = r.star_id

  UNION ALL
  -- Imported datasets (DACE, HARPS-RVBank, OHP, VizieR, uploads)
  SELECT sm.target, sm.matched_alias, sm.star_id,
         'Jump import: ' || concat_ws(' ', d.source, NULLIF(d.release, ''), NULLIF(d.variant, '')),
         d.slug || ':' || COALESCE(e.archive_row_id, e.id::text),
         concat_ws(' ', COALESCE(ei.code, ei.label), NULLIF(e.mode, '')),
         e.epoch_bjd, COALESCE(d.time_scale, 'unknown'), e.rv, e.rv_err,
         concat_ws('; ', 'drp_version=' || e.drp_version, 'rv_method=' || d.rv_method,
                   'ccf_mask=' || e.ccf_mask),
         d.corrections::text,
         concat_ws('; ', 'isjunk=' || e.isjunk, 'qc_flag=' || e.qc_flag, 'qc_note=' || e.qc_note,
                   'visible=' || e.visible, 'dataset_status=' || ds.status,
                   'preferred=' || ds.preferred, 'superseded_by=' || ds.superseded_by_id,
                   'retired_at=' || ds.retired_at),
         concat_ws('; ', 'snr=' || e.snr, 'texp=' || e.texp),
         e.s_index, NULL::bigint, e.program_code,
         -- DACE's published datasets keep the bibcode in `release`
         COALESCE(NULLIF(ds.paper_bibcode, ''), NULLIF(d.bibcode, ''),
                  CASE WHEN d.release ~ '^[0-9]{4}[A-Za-z]' THEN d.release END),
         'access=' || d.access,
         NULL, NULL::float8, d.rv_unit_published, e.texp
  FROM core_externalrv e
  JOIN core_rvdatasetstar ds ON ds.id = e.dataset_star_id
  JOIN core_rvdataset d ON d.id = ds.dataset_id
  LEFT JOIN core_externalinstrument ei ON ei.id = d.instrument_id
  JOIN star_map sm ON sm.any_star_id = e.star_id
)
SELECT l.target, l.matched_alias, s.name AS jump_star_name, l.database, l.obs_id,
       l.instrument, l.time_bjd, l.time_scale, l.rv_ms, l.rv_err_ms, l.reduction,
       l.corrections, l.quality_flag, l.counts_or_snr, l.s_index,
       COALESCE(NULLIF(concat_ws(' / ', p.progid, NULLIF(p.pi_name, '')), ''), l.program_code) AS program,
       l.publication, l.access, l.utctime, l.jd_kbcvel, l.rv_unit, l.exptime_s
FROM legs l
JOIN core_star s ON s.id = l.star_id
LEFT JOIN core_program p ON p.id = l.program_id
ORDER BY l.target, l.instrument, l.time_bjd;
