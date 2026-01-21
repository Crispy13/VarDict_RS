/// Debug test to compare CigarParser behavior with Java for specific positions
use rust_htslib::bam::{Read as BamRead, IndexedReader};
use vardict_rs::data::reference::Reference;
use vardict_rs::data::region::Region;
use vardict_rs::mods::cigar_parser::CigarParser;
use vardict_rs::scopedata::global_read_only_scope::{GlobalReadOnlyScope, INSTANCE};
use vardict_rs::conf::Configuration;
use std::sync::Arc;
use std::path::Path;

#[test]
#[ignore]
fn test_cigar_parser_position_168714() {
    // Setup logging
    use crackle_kit::tracing::level_filters::LevelFilter;
    crackle_kit::tracing_kit::setup_logging_stderr_only_verbose(LevelFilter::DEBUG);
    
    // Load reference for the region
    let ref_path = "/home/eck/workspace/VarDictJava/tests/integration/reference/hs37d5.fa";
    if !Path::new(ref_path).exists() {
        eprintln!("Reference file not found: {}", ref_path);
        return;
    }
    
    // Read reference sequence using rust-htslib
    // Note: fetch_seq_string uses 0-based coordinates, but our BED uses 1-based
    // BED region: 168527-168759 (1-based, inclusive)
    // For fetch: use 168526-168759 (0-based, half-open interval [start, end))
    let ref_reader = rust_htslib::faidx::Reader::from_path(ref_path)
        .expect("Failed to open reference");
    let ref_seq = ref_reader.fetch_seq_string("20", 168526, 168759)
        .expect("Failed to fetch reference sequence");
    
    // Reference stores 1-based start position
    let reference = Reference::from_seq_with_start(ref_seq.as_bytes(), 168527);
    
    println!("Reference loaded: {} bases", reference.ref_seq.len());
    println!("Reference at 168714: {:?}", reference.get(168714).map(|b| b as char));
    
    // Open BAM file
    let bam_path = "VarDictJava/tests/integration/input/NA12878.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam";
    if !Path::new(bam_path).exists() {
        eprintln!("BAM file not found: {}", bam_path);
        return;
    }
    
    let mut bam = IndexedReader::from_path(bam_path).expect("Failed to open BAM");
    
    // Create region
    let region = Region::new("20".to_string(), 168527, 168759, "test".to_string());
    
    // Fetch reads covering position 168714
    // Note: position 168714 is 1-based, BAM uses 0-based coordinates
    // So we need to fetch around 168713 in BAM coordinates
    bam.fetch(("20", 168600, 168800)).expect("Failed to fetch region");
    
    let mut reads_covering_168714 = Vec::new();
    let mut total_reads = 0;
    
    let mut t_reads_fwd = 0;
    let mut t_reads_rev = 0;
    let mut t_reads_fwd_oriented = 0;
    let mut t_reads_rev_oriented = 0;
    let mut t_reads_aligned: Vec<(String, u16, u8, String, char, u8, bool)> = Vec::new();
    let mut t_reads_oriented: Vec<(String, u16, u8, String, char, char, u8, bool)> = Vec::new();

    let complement = |base: char| -> char {
        match base {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        }
    };

    for result in bam.records() {
        total_reads += 1;
        let record = result.expect("Failed to read record");
        
        // Check if read covers position 168714 (1-based)
        // BAM pos() is 0-based, so record.pos() == 168713 means starts at 1-based 168714
        let start_0based = record.pos() as i64;
        let start_1based = start_0based + 1;
        let cigar = record.cigar();
        let end_0based = cigar.end_pos() as i64;  // exclusive
        let end_1based = end_0based;  // Because end is exclusive in 0-based, it's the same as 1-based inclusive-1
        
        // Check if position 168714 (1-based) is covered
        // In 0-based: starts <= 168713 && end > 168713
        if start_0based <= 168713 && end_0based > 168713 {
            // Get base at position 168714 (1-based) = 168713 (0-based)
            let ref_pos_in_read = (168713 - start_0based) as u32;
            
            // Try to get the read position from CIGAR
            // For simple matches, we can calculate directly
            let seq = record.seq();
            let qual = record.qual();
            
            // For a simple match CIGAR (e.g., 76M), the read position equals ref position
            // But for more complex CIGARs with indels, we need read_pos
            let read_pos_opt = cigar.read_pos(ref_pos_in_read, false, false);
            
            // Compute aligned base at 168714 using CIGAR walk
            let mut ref_pos = start_0based as i64;
            let mut read_pos = 0usize;
            let mut aligned_base: Option<(char, u8)> = None;
            for op in cigar.iter() {
                match *op {
                    rust_htslib::bam::record::Cigar::Match(l)
                    | rust_htslib::bam::record::Cigar::Equal(l)
                    | rust_htslib::bam::record::Cigar::Diff(l) => {
                        let l = l as i64;
                        if ref_pos <= 168713 && 168713 < ref_pos + l {
                            let offset = (168713 - ref_pos) as usize;
                            let idx = read_pos + offset;
                            if idx < seq.len() {
                                aligned_base = Some((seq.as_bytes()[idx] as char, qual[idx]));
                            }
                            break;
                        }
                        ref_pos += l;
                        read_pos += l as usize;
                    }
                    rust_htslib::bam::record::Cigar::Ins(l) => {
                        read_pos += l as usize;
                    }
                    rust_htslib::bam::record::Cigar::Del(l)
                    | rust_htslib::bam::record::Cigar::RefSkip(l) => {
                        ref_pos += l as i64;
                    }
                    rust_htslib::bam::record::Cigar::SoftClip(l) => {
                        read_pos += l as usize;
                    }
                    rust_htslib::bam::record::Cigar::HardClip(_) | rust_htslib::bam::record::Cigar::Pad(_) => {}
                }
            }

            let base_info = aligned_base;
            
            println!("\n=== Read {} (matches position 168714) ===", reads_covering_168714.len() + 1);
            println!("Name: {}", String::from_utf8_lossy(record.qname()));
            println!("Position: {}-{} (0-based: {}-{})", start_1based, end_1based, start_0based, end_0based);
            println!("Ref pos in read (0-based offset): {}", ref_pos_in_read);
            println!("Read pos result: {:?}", read_pos_opt);
            println!("Read length: {}", seq.len());
            
            println!("Strand: {}", if record.is_reverse() { "Reverse" } else { "Forward" });
            println!("CIGAR: {}", cigar.to_string());
            println!("MAPQ: {}", record.mapq());
            println!("Flags: 0x{:x}", record.flags());
            println!("Is_duplicate: {}", record.is_duplicate());
            println!("Is_secondary: {}", record.is_secondary());
            println!("Is_supplementary: {}", record.is_supplementary());
            
            if let Some((base, qual)) = base_info {
                println!("Base at 168714: {} (qual={})", base, qual);
                if base == 'T' {
                    if record.is_reverse() {
                        t_reads_rev += 1;
                    } else {
                        t_reads_fwd += 1;
                    }
                    t_reads_aligned.push((
                        String::from_utf8_lossy(record.qname()).to_string(),
                        record.flags(),
                        record.mapq(),
                        cigar.to_string(),
                        base,
                        qual,
                        record.is_reverse(),
                    ));
                }

                let oriented_base = if record.is_reverse() {
                    complement(base)
                } else {
                    base
                };
                if oriented_base == 'T' {
                    if record.is_reverse() {
                        t_reads_rev_oriented += 1;
                    } else {
                        t_reads_fwd_oriented += 1;
                    }
                    t_reads_oriented.push((
                        String::from_utf8_lossy(record.qname()).to_string(),
                        record.flags(),
                        record.mapq(),
                        cigar.to_string(),
                        base,
                        oriented_base,
                        qual,
                        record.is_reverse(),
                    ));
                }
            } else {
                println!("Base at 168714: N/A (not in aligned region)");
            }
            
            // Try to get NM tag
            if let Ok(nm) = record.aux(b"NM") {
                println!("NM tag: {:?}", nm);
            }
            
            reads_covering_168714.push(record);
        }
    }
    
    println!("\n=== Summary ===");
    println!("Total reads fetched: {}", total_reads);
    println!("Reads covering position 168714: {}", reads_covering_168714.len());
    println!("Reads with T at 168714 (CIGAR-aligned): fwd={}, rev={}", t_reads_fwd, t_reads_rev);
    println!("Reads with T at 168714 (oriented to reference): fwd={}, rev={}", t_reads_fwd_oriented, t_reads_rev_oriented);

    println!("\n=== Reads with T at 168714 (CIGAR-aligned) ===");
    for (qname, flags, mapq, cigar, base, qual, is_reverse) in &t_reads_aligned {
        println!(
            "{}\tflags=0x{:x}\tmapq={}\tcigar={}\tbase={}\tqual={}\trev={}",
            qname, flags, mapq, cigar, base, qual, is_reverse
        );
    }

    println!("\n=== Reads with T at 168714 (oriented to reference) ===");
    for (qname, flags, mapq, cigar, base, oriented_base, qual, is_reverse) in &t_reads_oriented {
        println!(
            "{}\tflags=0x{:x}\tmapq={}\tcigar={}\tbase={}\toriented={}\tqual={}\trev={}",
            qname, flags, mapq, cigar, base, oriented_base, qual, is_reverse
        );
    }
    
    // Now parse these reads with CigarParser
    println!("\n=== Parsing reads with CigarParser ===");
    
    // Create and initialize global scope
    let mut conf = Configuration::default();
    conf.goodq = 22.5;
    conf.vext = 2;
    conf.disable_sv = true;
    conf.perform_local_realignment = true;
    let mut scope = GlobalReadOnlyScope::default();
    scope.conf = conf;
    // Add chr_lens for chromosome 20 (hs37d5 chromosome 20 length is 63025520)
    scope.chr_lens.insert("20".to_string(), 63025520);
    
    // Initialize the global INSTANCE
    let _ = INSTANCE.set(scope.clone());
    let instance = Arc::new(scope);
    
    let mut parser = CigarParser::new(region.clone(), reference.clone(), instance);
    
    // Parse the records
    let mut records_for_parsing = reads_covering_168714;
    parser.process_records(records_for_parsing.iter_mut())
        .expect("Failed to parse records");
    
    println!("Parsed {} reads successfully", records_for_parsing.len());
    
    // Check what variants were created at position 168714
    println!("\n=== Variants at position 168714 ===");
    
    // Check both position 168714 directly and nearby positions
    for pos in 168710..168720 {
        if let Some(var_map) = parser.get_non_insertion_vars().get(&pos) {
            println!("\nVariants at position {}:", pos);
            for (desc, variant) in var_map.iter() {
                println!("  {:?}: alt_depth={} (fwd={}, rev={})",
                    desc,
                    variant.alt_depth,
                    variant.alt_depth_fwd,
                    variant.alt_depth_rev);
            }
        }
    }
    
    // Compare with expected Java output:
    // Java: alt_depth=2 (fwd=1, rev=1)
    // We should see the same
    println!("\n=== Expected (from Java) ===");
    println!("Position 168714, C→T: alt_depth=2 (fwd=1, rev=1)");
    println!("\nIf Rust shows alt_depth=1 with only forward or only reverse,");
    println!("then one read is not being counted properly.");

}

#[test]
#[ignore]
fn test_cigar_parser_variant_76749_matches_java() {
    use crackle_kit::tracing::level_filters::LevelFilter;
    use rust_htslib::bam::Record;
    use vardict_rs::data::bam_reader::BamReader;
    use vardict_rs::mods::vardict_pipeline::VarDictPipeline;
    use vardict_rs::variants::variants::VarDesc;

    crackle_kit::tracing_kit::setup_logging_stderr_only_verbose(LevelFilter::DEBUG);

    let ref_path = "/home/eck/workspace/VarDictJava/tests/integration/reference/hs37d5.fa";
    let bam_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/input/NA12878.chrom20.ILLUMINA.bwa.CEU.exome.20121211.bam";
    let java_output_path = "/home/eck/workspace/vardict_rs/VarDictJava/tests/integration/raw_input/raw.vardict.simple.chr20.nosv.var";
    if !Path::new(ref_path).exists() || !Path::new(bam_path).exists() {
        eprintln!("Missing reference or BAM file");
        return;
    }
    if !Path::new(java_output_path).exists() {
        eprintln!("Missing Java output file: {}", java_output_path);
        return;
    }

    // Region around the first mismatch: 20:76646-76845 (1-based, inclusive)
    let region = Region::new("20".to_string(), 76646, 76845, "test".to_string());

    // Load reference sequence (0-based, half-open)
    let ref_reader = rust_htslib::faidx::Reader::from_path(ref_path)
        .expect("Failed to open reference");
    let ref_seq = ref_reader
        .fetch_seq_string("20", 76645, 76845)
        .expect("Failed to fetch reference sequence");
    let reference = Reference::from_seq_with_start(ref_seq.as_bytes(), 76646);

    // Configure global scope to match simple mode defaults
    let mut conf = Configuration::default();
    conf.goodq = 22.5;
    conf.vext = 2;
    conf.disable_sv = true;
    conf.perform_local_realignment = true;

    let mut scope = GlobalReadOnlyScope::default();
    scope.conf = conf;
    scope.chr_lens.insert("20".to_string(), 63025520);
    let _ = INSTANCE.set(scope.clone());
    let instance = Arc::new(scope);

    let pipeline = VarDictPipeline::new("abc");
    let sam_filter = instance.conf.sam_filter;
    let mut bam_reader = BamReader::open(bam_path).expect("Failed to open BAM");
    bam_reader
        .fetch(region.chr(), region.start(), region.end())
        .expect("Failed to fetch BAM region");

    let mut records = Vec::new();
    let mut record = Record::new();
    while bam_reader.read(&mut record).unwrap_or(false) {
        // Mirror VarDictPipeline::passes_preprocess
        if sam_filter != 0 && (record.flags() & (sam_filter as u16)) != 0 {
            continue;
        }
        if pipeline.min_mapping_quality > 0 && record.mapq() < pipeline.min_mapping_quality {
            continue;
        }
        const SECONDARY_ALIGNMENT: u16 = 0x100;
        if (record.flags() & SECONDARY_ALIGNMENT) != 0 && sam_filter != 0 {
            continue;
        }
        let seq = record.seq();
        if seq.len() == 0 || (seq.len() == 1 && seq.as_bytes()[0] == b'*') {
            continue;
        }
        records.push(record.clone());
    }

    let mut parser = CigarParser::new(region.clone(), reference.clone(), instance);
    parser
        .process_records(records.iter_mut())
        .expect("Failed to parse records");

    #[derive(Debug)]
    struct JavaVariantCounts {
        pos: i64,
        ref_base: u8,
        alt_base: u8,
        alt_depth: usize,
        alt_fwd: usize,
        alt_rev: usize,
        ref_fwd: usize,
        ref_rev: usize,
        hicov: usize,
    }

    fn extract_java_variant(java_output_path: &str) -> Option<JavaVariantCounts> {
        let content = std::fs::read_to_string(java_output_path).ok()?;
        for line in content.lines() {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 30 {
                continue;
            }
            let chr = parts.get(2)?;
            let start = parts.get(3)?;
            let end = parts.get(4)?;
            let ref_base = parts.get(5)?;
            let alt_base = parts.get(6)?;
            if *chr == "20" && *start == "76749" && *end == "76749" && *ref_base == "A" && *alt_base == "G" {
                let alt_depth = parts.get(8)?.parse().ok()?;
                let ref_fwd = parts.get(9)?.parse().ok()?;
                let ref_rev = parts.get(10)?.parse().ok()?;
                let alt_fwd = parts.get(11)?.parse().ok()?;
                let alt_rev = parts.get(12)?.parse().ok()?;
                let hicov = parts.get(29)?.parse().ok()?;

                return Some(JavaVariantCounts {
                    pos: 76749,
                    ref_base: ref_base.as_bytes()[0],
                    alt_base: alt_base.as_bytes()[0],
                    alt_depth,
                    alt_fwd,
                    alt_rev,
                    ref_fwd,
                    ref_rev,
                    hicov,
                });
            }
        }
        None
    }

    let java_counts = extract_java_variant(java_output_path)
        .expect("Failed to extract Java variant counts for 20:76749 A>G");

    println!("Java extracted counts: {:?}", java_counts);

    let pos = java_counts.pos;
    let ref_base = reference.get(pos).unwrap_or(b'N');
    let alt_base = java_counts.alt_base;

    let var_map = parser
        .get_non_insertion_vars()
        .get(&pos)
        .expect("No variants found at 76749");

    let ref_var = var_map
        .get(&VarDesc::SNV { ref_base })
        .expect("Missing reference variant at 76749");
    let alt_key = VarDesc::SNV { ref_base: alt_base };
    let alt_var = match var_map.get(&alt_key) {
        Some(v) => v,
        None => {
            println!("Alt variant missing at 76749. Available keys:");
            for (desc, variant) in var_map.iter() {
                println!("  {:?}: alt_depth={} fwd={} rev={}", desc, variant.alt_depth, variant.alt_depth_fwd, variant.alt_depth_rev);
            }
            println!("Nearby positions with alt base '{}' in non_insertion_vars:", alt_base as char);
            let start_pos = pos.saturating_sub(5);
            let end_pos = pos + 5;
            for p in start_pos..=end_pos {
                if let Some(nearby_map) = parser.get_non_insertion_vars().get(&p) {
                    let key = VarDesc::SNV { ref_base: alt_base };
                    if let Some(variant) = nearby_map.get(&key) {
                        println!(
                            "  pos {}: alt_depth={} fwd={} rev={}",
                            p,
                            variant.alt_depth,
                            variant.alt_depth_fwd,
                            variant.alt_depth_rev
                        );
                    }
                }
            }
            panic!("Missing alt variant at 76749");
        }
    };

    println!("Ref {} counts: alt_depth={} fwd={} rev={}",
        ref_base as char,
        ref_var.alt_depth,
        ref_var.alt_depth_fwd,
        ref_var.alt_depth_rev,
    );
    println!("Alt {} counts: alt_depth={} fwd={} rev={}",
        alt_base as char,
        alt_var.alt_depth,
        alt_var.alt_depth_fwd,
        alt_var.alt_depth_rev,
    );

    assert_eq!(alt_var.alt_depth, java_counts.alt_depth, "alt_depth mismatch for 76749");
    assert_eq!(alt_var.alt_depth_fwd, java_counts.alt_fwd, "alt fwd mismatch for 76749");
    assert_eq!(alt_var.alt_depth_rev, java_counts.alt_rev, "alt rev mismatch for 76749");

    assert_eq!(ref_var.alt_depth_fwd, java_counts.ref_fwd, "ref fwd mismatch for 76749");
    assert_eq!(ref_var.alt_depth_rev, java_counts.ref_rev, "ref rev mismatch for 76749");
}
