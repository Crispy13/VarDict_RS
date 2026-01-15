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
