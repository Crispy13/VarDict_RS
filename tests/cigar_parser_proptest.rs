use std::sync::{Arc, LazyLock, Mutex, MutexGuard};

use proptest::prelude::*;
use proptest::test_runner::Config as ProptestConfig;
use rust_htslib::bam::{
    Record,
    record::{Cigar, CigarString},
};
use vardict_rs::conf::Configuration;
use vardict_rs::data::reference::Reference;
use vardict_rs::data::region::Region;
use vardict_rs::mods::cigar_parser::CigarParser;
use vardict_rs::prelude::SmallVecBytes;
use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
use vardict_rs::variants::variants::VarDesc;

static TEST_SCOPE_LOCK: LazyLock<Mutex<()>> = LazyLock::new(|| Mutex::new(()));

#[derive(Debug, Clone)]
struct ParserCase {
    ops: Vec<(char, u32)>,
    query_seq: Vec<u8>,
    is_reverse: bool,
}

fn lock_test_scope() -> MutexGuard<'static, ()> {
    TEST_SCOPE_LOCK.lock().expect("test scope mutex poisoned")
}

#[allow(invalid_reference_casting)]
fn install_test_scope(scope: GlobalReadOnlyScope) {
    if let Some(existing) = INSTANCE.get() {
        unsafe {
            let ptr = existing as *const GlobalReadOnlyScope as *mut GlobalReadOnlyScope;
            let _ = std::ptr::replace(ptr, scope);
        }
    } else {
        let _ = INSTANCE.set(scope);
    }
}

fn parser_scope() -> GlobalReadOnlyScope {
    let mut conf = Configuration::default();
    conf.disable_sv = true;
    conf.perform_local_realignment = false;
    conf.sam_filter = 0;

    let mut scope = GlobalReadOnlyScope::default();
    scope.conf = conf;
    scope
}

fn arb_base() -> impl Strategy<Value = u8> {
    prop_oneof![Just(b'A'), Just(b'C'), Just(b'G'), Just(b'T')]
}

fn arb_dna_seq(len: usize) -> impl Strategy<Value = Vec<u8>> {
    prop::collection::vec(arb_base(), len)
}

fn arb_general_cigar_op() -> impl Strategy<Value = (char, u32)> {
    (
        prop_oneof![
            Just('M'),
            Just('I'),
            Just('D'),
            Just('S'),
            Just('H'),
            Just('N'),
            Just('='),
            Just('X')
        ],
        1u32..=64,
    )
}

fn arb_cigar_ops() -> impl Strategy<Value = Vec<(char, u32)>> {
    prop::collection::vec(arb_general_cigar_op(), 1..=8).prop_filter(
        "cigar must consume at least one query or reference base",
        |ops| cigar_query_consumed(ops) > 0 || cigar_ref_consumed(ops) > 0,
    )
}

fn arb_alignment_op() -> impl Strategy<Value = (char, u32)> {
    (prop_oneof![Just('M'), Just('='), Just('X')], 1u32..=12)
}

fn arb_body_op() -> impl Strategy<Value = (char, u32)> {
    (
        prop_oneof![Just('M'), Just('I'), Just('D'), Just('='), Just('X')],
        1u32..=12,
    )
}

fn arb_clip_op() -> impl Strategy<Value = Option<(char, u32)>> {
    prop_oneof![
        Just(None),
        (prop_oneof![Just('S'), Just('H')], 1u32..=8).prop_map(Some),
    ]
}

fn arb_parser_case() -> impl Strategy<Value = ParserCase> {
    (
        arb_clip_op(),
        arb_alignment_op(),
        prop::collection::vec(arb_body_op(), 0..=3),
        arb_alignment_op(),
        arb_clip_op(),
        any::<bool>(),
    )
        .prop_flat_map(|(lead, first, middle, last, tail, is_reverse)| {
            let mut ops = Vec::with_capacity(middle.len() + 4);
            if let Some(lead) = lead {
                ops.push(lead);
            }
            ops.push(first);
            ops.extend(middle);
            ops.push(last);
            if let Some(tail) = tail {
                ops.push(tail);
            }

            let query_len = cigar_query_consumed(&ops);
            (Just(ops), arb_dna_seq(query_len), Just(is_reverse)).prop_map(
                |(ops, query_seq, is_reverse)| ParserCase {
                    ops,
                    query_seq,
                    is_reverse,
                },
            )
        })
}

fn arb_var_desc() -> impl Strategy<Value = VarDesc> {
    prop_oneof![
        arb_base().prop_map(VarDesc::snv_key),
        prop::collection::vec(arb_base(), 1..=16).prop_map(|seq| VarDesc::insertion(&seq)),
        (1u32..=256).prop_map(VarDesc::deletion),
        (
            prop::collection::vec(arb_base(), 1..=8),
            prop::collection::vec(arb_base(), 1..=8),
        )
            .prop_filter(
                "complex ref and alt must differ",
                |(ref_seq, alt_seq)| ref_seq != alt_seq
            )
            .prop_map(|(ref_seq, alt_seq)| VarDesc::complex(&ref_seq, &alt_seq)),
        prop::collection::vec(prop_oneof![Just(b'-'), arb_base()], 1..=16).prop_map(|desc| {
            VarDesc::Raw {
                desc: SmallVecBytes::from_slice(&desc),
            }
        }),
    ]
}

fn cigar_to_string(ops: &[(char, u32)]) -> String {
    ops.iter()
        .map(|(op, len)| format!("{}{}", len, op))
        .collect::<String>()
}

fn parse_cigar_string(cigar: &str) -> Vec<Cigar> {
    let mut len = 0u32;
    let mut ops = Vec::new();

    for byte in cigar.bytes() {
        if byte.is_ascii_digit() {
            len = len * 10 + u32::from(byte - b'0');
            continue;
        }

        assert!(len > 0, "operation length must be positive");
        ops.push(op_to_cigar(byte as char, len));
        len = 0;
    }

    assert_eq!(len, 0, "unterminated cigar length");
    ops
}

fn op_to_cigar(op: char, len: u32) -> Cigar {
    match op {
        'M' => Cigar::Match(len),
        'I' => Cigar::Ins(len),
        'D' => Cigar::Del(len),
        'S' => Cigar::SoftClip(len),
        'H' => Cigar::HardClip(len),
        'N' => Cigar::RefSkip(len),
        '=' => Cigar::Equal(len),
        'X' => Cigar::Diff(len),
        _ => panic!("unsupported CIGAR op: {op}"),
    }
}

fn cigar_query_consumed(ops: &[(char, u32)]) -> usize {
    ops.iter()
        .map(|(op, len)| match op {
            'M' | 'I' | 'S' | '=' | 'X' => *len as usize,
            'D' | 'H' | 'N' => 0,
            _ => 0,
        })
        .sum()
}

fn cigar_ref_consumed(ops: &[(char, u32)]) -> usize {
    ops.iter()
        .map(|(op, len)| match op {
            'M' | 'D' | 'N' | '=' | 'X' => *len as usize,
            'I' | 'S' | 'H' => 0,
            _ => 0,
        })
        .sum()
}

fn mismatch_base(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        _ => b'A',
    }
}

fn build_reference_from_query(ops: &[(char, u32)], query_seq: &[u8]) -> Vec<u8> {
    let mut query_index = 0usize;
    let mut reference = Vec::with_capacity(cigar_ref_consumed(ops));

    for (op, len) in ops {
        match op {
            'M' | '=' => {
                for _ in 0..*len as usize {
                    let base = query_seq[query_index];
                    reference.push(base);
                    query_index += 1;
                }
            }
            'X' => {
                for _ in 0..*len as usize {
                    let base = query_seq[query_index];
                    reference.push(mismatch_base(base));
                    query_index += 1;
                }
            }
            'I' | 'S' => {
                query_index += *len as usize;
            }
            'D' | 'N' => {
                reference.extend(std::iter::repeat_n(b'A', *len as usize));
            }
            'H' => {}
            _ => unreachable!("unsupported parser op"),
        }
    }

    reference
}

fn build_record(ops: &[(char, u32)], query_seq: &[u8], is_reverse: bool, pos: i64) -> Record {
    let mut record = Record::new();
    let cigar = CigarString(
        ops.iter()
            .map(|(op, len)| op_to_cigar(*op, *len))
            .collect::<Vec<_>>(),
    );
    let quals = vec![40u8; query_seq.len()];

    record.set(b"prop_read", Some(&cigar), query_seq, &quals);
    record.set_tid(0);
    record.set_pos(pos);
    record.set_mapq(60);
    record.set_flags(if is_reverse { 0x10 } else { 0 });
    record.set_mtid(0);
    record.set_mpos(-1);
    record.set_insert_size(0);

    record
}

fn run_parser_case(case: &ParserCase) -> anyhow::Result<(CigarParser, Region, i64, usize)> {
    let _guard = lock_test_scope();

    let scope = parser_scope();
    install_test_scope(scope.clone());

    let region_start = 1_000usize;
    let alignment_start = region_start as i64 + 5;
    let ref_consumed = cigar_ref_consumed(&case.ops);
    let region_end = alignment_start as usize + ref_consumed + 10;
    let region = Region::new(
        "prop_chr".to_string(),
        region_start,
        region_end,
        "PROP".to_string(),
    );

    let ref_segment = build_reference_from_query(&case.ops, &case.query_seq);
    let mut reference_seq = vec![b'A'; region.len()];
    let start_offset = alignment_start as usize - region_start;
    let end_offset = start_offset + ref_segment.len();
    if !ref_segment.is_empty() {
        reference_seq[start_offset..end_offset].copy_from_slice(&ref_segment);
    }

    let reference = Reference::new_with_start(reference_seq, region_start as i64);
    let mut parser = CigarParser::new(region.clone(), reference, Arc::new(scope));
    let mut record = build_record(
        &case.ops,
        &case.query_seq,
        case.is_reverse,
        alignment_start - 1,
    );

    parser.process_record(&mut record)?;

    Ok((parser, region, alignment_start, ref_consumed))
}

proptest! {
    #![proptest_config(ProptestConfig {
        cases: 1000,
        failure_persistence: None,
        .. ProptestConfig::default()
    })]

    #[test]
    fn prop_cigar_consumed_bases_match_operation_semantics(ops in arb_cigar_ops()) {
        let cigar_string = cigar_to_string(&ops);
        let parsed = parse_cigar_string(&cigar_string);

        let expected_query = cigar_query_consumed(&ops);
        let expected_ref = cigar_ref_consumed(&ops);
        let parsed_query: usize = parsed
            .iter()
            .map(|op| match op {
                Cigar::Match(len) | Cigar::Ins(len) | Cigar::SoftClip(len) | Cigar::Equal(len) | Cigar::Diff(len) => *len as usize,
                Cigar::Del(_) | Cigar::HardClip(_) | Cigar::RefSkip(_) => 0,
                _ => 0,
            })
            .sum();
        let parsed_ref: usize = parsed
            .iter()
            .map(|op| match op {
                Cigar::Match(len) | Cigar::Del(len) | Cigar::RefSkip(len) | Cigar::Equal(len) | Cigar::Diff(len) => *len as usize,
                Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => 0,
                _ => 0,
            })
            .sum();

        prop_assert_eq!(parsed_query, expected_query);
        prop_assert_eq!(parsed_ref, expected_ref);
    }

    #[test]
    fn prop_cigar_string_round_trips_through_rust_htslib(ops in arb_cigar_ops()) {
        let cigar_string = cigar_to_string(&ops);
        let parsed = parse_cigar_string(&cigar_string);
        let round_trip = CigarString(parsed).to_string();

        prop_assert_eq!(round_trip, cigar_string);
    }

    #[test]
    fn prop_reference_get_stays_in_bounds(
        seq in prop::collection::vec(arb_base(), 1..=256),
        region_start in 1i64..=10_000,
        probe_seed in any::<usize>(),
    ) {
        let reference = Reference::new_with_start(seq.clone(), region_start);
        let valid_index = probe_seed % seq.len();
        let valid_pos = region_start + valid_index as i64;

        prop_assert_eq!(reference.get(valid_pos), Some(seq[valid_index]));
        prop_assert_eq!(reference.get(region_start - 1), None);
        prop_assert_eq!(reference.get(region_start + seq.len() as i64), None);
    }

    #[test]
    fn prop_var_desc_key_and_classification_methods_stay_consistent(desc in arb_var_desc()) {
        let key = desc.to_key_string();

        prop_assert!(desc.key_equals(&key));
        prop_assert_eq!(desc.is_deletion(), key.starts_with('-'));

        match &desc {
            VarDesc::SNV { ref_base } => {
                prop_assert_eq!(desc.variant_type(), "SNV");
                prop_assert!(desc.key_matches_byte(*ref_base));
            }
            VarDesc::Ins { .. } => {
                prop_assert_eq!(desc.variant_type(), "Insertion");
                prop_assert!(key.starts_with('+'));
            }
            VarDesc::Del { .. } => {
                prop_assert_eq!(desc.variant_type(), "Deletion");
                prop_assert!(desc.is_deletion());
            }
            VarDesc::Complex { .. } => {
                prop_assert_eq!(desc.variant_type(), "Complex");
                prop_assert!(key.contains('>'));
            }
            VarDesc::Raw { .. } => {
                prop_assert_eq!(desc.variant_type(), "Raw");
            }
        }
    }

    #[test]
    fn prop_cigar_parser_accepts_well_formed_single_record_inputs(case in arb_parser_case()) {
        let result = run_parser_case(&case);
        prop_assert!(result.is_ok(), "parser rejected case: {:?}", case);
    }

    #[test]
    fn prop_cigar_parser_outputs_stay_within_expected_bounds(case in arb_parser_case()) {
        let (parser, region, alignment_start, ref_consumed) = run_parser_case(&case)
            .expect("parser case should succeed");
        let region_start = region.start() as i64;
        let region_end = region.end() as i64;
        let coverage_min = alignment_start - 1;
        let coverage_max = alignment_start
            + cigar_query_consumed(&case.ops).max(ref_consumed) as i64;

        for pos in parser.get_non_insertion_vars().keys() {
            prop_assert!((*pos >= region_start) && (*pos <= region_end));
        }

        for pos in parser.get_insertion_vars().keys() {
            prop_assert!((*pos >= region_start) && (*pos <= region_end));
        }

        for (pos, cov) in parser.get_ref_coverage().iter_sorted() {
            prop_assert!(cov > 0);
            prop_assert!((pos >= coverage_min) && (pos <= coverage_max));
        }
    }
}
