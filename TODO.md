## [x] Make test functions using extra test data ignored
- Added `#[ignore = "requires VarDictJava integration test data (not in repo)"]` to:
  1. `src/mods/cigar_parser.rs` — `test_cigar_parser_non_insertion_variants_first_bed_region`
  2. `src/mods/vardict_pipeline.rs` — `test_record_preprocessor_dump_all_bed_regions`
  3. `tests/variant_realigner_fixture_test.rs` — `test_variant_realigner_snapshot_fixtures`
  4. `tests/tovars_fixture_test.rs` — `test_tovars_builder_snapshot_fixture`

## [ ] Paired mode bam
1. https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/superseded-2022-data/BCM_ILMN-somatic-analysis_20220816/ (WGS, ~150GB)
2. https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES/
- Tumor: WES_IL_T_1.bwa.dedup.bam (Illumina replicate 1)
- Normal: WES_IL_N_1.bwa.dedup.bam (Illumina replicate 1)
3. https://ftp.ncbi.nlm.nih.gov/?