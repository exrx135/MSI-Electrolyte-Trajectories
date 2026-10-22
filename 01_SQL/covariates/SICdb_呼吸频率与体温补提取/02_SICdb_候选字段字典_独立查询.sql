-- Standalone PostgreSQL query. Native SICdb tables only; no earlier study tables required.
-- Only frozen 3831 cases; [sepsis onset, onset + 86400 seconds).
-- Read-only native sources; all generated objects are temporary.
CREATE TEMP TABLE IF NOT EXISTS sic_vital_session(dummy integer);
CREATE OR REPLACE FUNCTION pg_temp.t1_sic_native(out_table text, native_table text,
 fields jsonb, required boolean DEFAULT true, use_window boolean DEFAULT false)
RETURNS void LANGUAGE plpgsql AS $native$
DECLARE ns text; tbl text; candidates integer; f jsonb; col text; expr text := '';
 case_col text; time_col text; data_col text; sql_text text;
BEGIN
 SELECT count(DISTINCT n.nspname),min(n.nspname) INTO candidates,ns
 FROM pg_namespace n WHERE n.nspname NOT LIKE 'pg_%' AND n.nspname <> 'information_schema'
 AND EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND lower(c.relname)='cases' AND c.relkind IN ('r','p','v','m','f'))
 AND EXISTS(SELECT 1 FROM pg_class c WHERE c.relnamespace=n.oid AND lower(c.relname)='d_references' AND c.relkind IN ('r','p','v','m','f'));
 IF candidates<>1 THEN RAISE EXCEPTION 'Expected one native SICdb schema containing cases+d_references; found %',candidates; END IF;
 SELECT c.relname INTO tbl FROM pg_class c JOIN pg_namespace n ON n.oid=c.relnamespace
 WHERE n.nspname=ns AND lower(c.relname)=native_table AND c.relkind IN ('r','p','v','m','f');
 IF tbl IS NULL THEN
  IF required THEN RAISE EXCEPTION 'Native SICdb source %.% missing',ns,native_table; END IF;
  RETURN;
 END IF;
 FOR f IN SELECT value FROM jsonb_array_elements(fields) LOOP
  SELECT a.attname INTO col FROM pg_attribute a
  WHERE a.attrelid=format('%I.%I',ns,tbl)::regclass AND a.attnum>0 AND NOT a.attisdropped
  AND regexp_replace(lower(a.attname),'[^a-z0-9]','','g') IN (SELECT jsonb_array_elements_text(f->'aliases'))
  ORDER BY (SELECT min(z.ordinality) FROM jsonb_array_elements_text(f->'aliases') WITH ORDINALITY z(value,ordinality)
    WHERE z.value=regexp_replace(lower(a.attname),'[^a-z0-9]','','g')) LIMIT 1;
  IF col IS NULL AND COALESCE((f->>'required')::boolean,false) THEN
   RAISE EXCEPTION 'Native %.% missing required field %',ns,tbl,f->>'name'; END IF;
  IF expr<>'' THEN expr:=expr||','; END IF;
  expr:=expr||CASE WHEN col IS NULL THEN format('NULL::%s AS %I',f->>'type',f->>'name')
    ELSE format('d.%I::%s AS %I',col,f->>'type',f->>'name') END;
  IF f->>'name'='case_id' THEN case_col:=col; END IF;
  IF f->>'name'='offset_sec' THEN time_col:=col; END IF;
 IF f->>'name'='data_id' THEN data_col:=col; END IF;
 END LOOP;
 IF native_table='cases' THEN expr:=expr||',to_jsonb(d) AS raw_record'; END IF;
 EXECUTE format('DROP TABLE IF EXISTS pg_temp.%I',out_table);
 sql_text:=format('CREATE TEMP TABLE %I AS SELECT %s FROM %I.%I d',out_table,expr,ns,tbl);
 IF case_col IS NOT NULL THEN sql_text:=sql_text||format(' JOIN t1_seed s ON d.%I::bigint=s.unit_id',case_col); END IF;
 IF native_table='data_float_h' THEN sql_text:=sql_text||format(' JOIN sic_signal_candidates c ON d.%I=c.ref_id',data_col); END IF;
 IF use_window THEN
  IF case_col IS NULL OR time_col IS NULL THEN RAISE EXCEPTION 'Time window columns not resolved'; END IF;
  sql_text:=sql_text||format(' WHERE d.%I::numeric>=s.onset_offset AND d.%I::numeric<s.onset_offset+86400',time_col,time_col);
 END IF;
 EXECUTE sql_text;
 IF case_col IS NOT NULL THEN EXECUTE format('CREATE INDEX ON %I(case_id)',out_table); END IF;
 RAISE NOTICE 'Native source: %.% -> current session %',ns,tbl,out_table;
END $native$;

DO $load$ BEGIN PERFORM pg_temp.t1_sic_native('sic_src_refs','d_references','[{"name": "ref_id", "aliases": ["referenceglobalid"], "type": "integer", "required": true}, {"name": "ref_value", "aliases": ["referencevalue"], "type": "text", "required": true}, {"name": "ref_name", "aliases": ["referencename"], "type": "text", "required": false}, {"name": "ref_unit", "aliases": ["referenceunit"], "type": "text", "required": false}]'::jsonb,true,false); END $load$;
DROP TABLE IF EXISTS pg_temp.sic_signal_candidates;
CREATE TEMP TABLE sic_signal_candidates AS
SELECT *, CASE WHEN concat_ws(' ',ref_value,ref_name) ~* '(temp|temperatur)' THEN 'temperature' ELSE 'respiratory' END AS candidate_family
FROM sic_src_refs
WHERE concat_ws(' ',ref_value,ref_name) ~* '(resp|breath|atem|frequen|temp|ventilat)' OR ref_id=709;
-- Broad matching intentionally retains measured, spontaneous, total and SET channels separately.

SELECT ref_id AS source_id,candidate_family,ref_value AS source_label,ref_name AS source_group,ref_unit AS native_unit,
candidate_family||'_'||ref_id||'_mean' AS mean_column,
candidate_family||'_'||ref_id||'_n' AS count_column
FROM sic_signal_candidates ORDER BY candidate_family,ref_id;
