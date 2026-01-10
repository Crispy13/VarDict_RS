# VarDictJava Code Examination - Summary Report

## Examination Complete ✓

I've thoroughly explored the VarDictJava codebase to understand the logic and architecture for porting to Rust.

---

## Key Findings

### 1. **Architecture Overview**
VarDict is a **pipeline-based variant caller** that processes genomic regions sequentially:
- Reads from BAM files
- Parses CIGAR strings to find variants
- Performs soft-clip realignment (novel feature)
- Calculates comprehensive statistics for each variant
- Applies quality filters
- Outputs in VCF/text format

### 2. **Three Main Processing Paths**

#### Path A: **Data Parsing** (Modules)
```
BAM → SAMFileParser → RecordPreprocessor → CigarParser → VariationRealigner
                                                            ↓
                                                      Low-level variations
```

#### Path B: **Variant Creation** (ToVarsBuilder) ⭐⭐⭐
```
Variations → Group by position → Calculate statistics → Create Variant objects
                                  ├─ Allele frequency
                                  ├─ Position in read
                                  ├─ Quality distribution
                                  ├─ Mapping quality
                                  ├─ Strand bias
                                  └─ MSI detection
```

#### Path C: **Output** (SimplePostProcessModule)
```
Variants → isGoodVar() filter → Adjust complex → Output
           ├─ Not strong bias
           ├─ At 2+ positions in reads
           ├─ At 2+ quality levels
           └─ Passes other thresholds
```

### 3. **Core Data Structures**

#### Variant.java (412 lines)
The **most important** class - holds all variant information:
- Identity: `refallele`, `varallele`, `vartype`
- Counts: `varsCountOnForward/Reverse`, `frequency`
- Quality metrics: `meanQuality`, `meanPosition`, `meanMappingQuality`
- Flags: `strandBiasFlag`, `isAtLeastAt2Positions`, `hasAtLeast2DiffQualities`
- Context: `leftseq` (20bp), `rightseq` (20bp)
- Specialized: `msi`, `msint`, `shift3`, `crispr`

#### Scope<T> (Generic Pipeline Container)
Carries data through pipeline stages:
```
Scope<InitialData> → Scope<VariationData> → Scope<RealignedVariationData> 
                  → Scope<AlignedVarsData>
```

#### Vars.java
Collection at each position:
- `variants: List<Variant>` - the variants found
- `referenceVariant: Variant` - reference call (for pileup)
- `sv: StructuralVariantFlags` - SV info

### 4. **Key Algorithms**

#### Soft-Clip Realignment (Novel Feature) ⭐
```
Read: ACGTACGT     (reference matches at 10-13, then 4 bases soft-clipped)
CIGAR: 10M4S

Actually:          ACGTACGTAAAA
Ref at region:  ACGTACGT----GGGG
                                ↑
                            deletion!
```
Without realignment: Soft clip not flagged as variant
With realignment: Deletion discovered in clipped sequence

#### Statistics Calculation
For each variant, calculate:
1. **Positions in reads** - List all positions where variant appears across reads
2. **Qualities in reads** - List all base qualities
3. **Strand counts** - Forward vs reverse
4. **Check diversity** - Is it at 2+ positions? 2+ quality levels?
5. **Calculate aggregates** - Mean, standard deviation

#### Quality Filter (isGoodVar)
```
PASS if:
✓ Not strong strand bias (unless in special region)
✓ Variant at 2+ positions across reads (prevents single-read errors)
✓ Variant at 2+ quality levels (prevents uniform-quality errors)
✓ Frequency ≥ minFreq
✓ Count ≥ minAltCount
✓ Not in MSI region OR MSI thresholds met
```

### 5. **Processing Pipeline (Simple Mode)**

**Sequential for each region:**
```
Region → BAM lookup → Parse reads → Parse CIGAR → Realign clips → 
Create variants → Calculate stats → Filter quality → Output
```

**Parallelizable for multiple regions:**
```
Region 1 ─┐
Region 2 ─┼─→ Thread Pool → Process independently → Queue output
Region 3 ─┤
Region 4 ─┘
```

### 6. **Output Format**

**Header:**
```
Sample Gene Chr Start End Ref Alt Depth AltDepth RefFwdReads RefRevReads 
AltFwdReads AltRevReads Genotype AF Bias PMean PStd QMean QStd MQ 
Sig_Noise HiAF ExtraAF shift3 MSI MSI_NT NM HiCnt HiCov 5pFlankSeq 
3pFlankSeq Seg VarType Duprate SV_info [CRISPR]
```

**Example data row:**
```
sample1 EGFR chr7 55086707 55086707 G A 150 45 20 25 22 23 HET 0.30 0 
45.5 12.3 28.4 3.2 55 SNP_NOISE 0.32 0.01 0 1.2 2 0 0.3 40 120 
ACGTACGTACGTACGTACGT ACGTACGTACGTACGTACGT chr7:55000000-55200000 SNP 
0.02 .
```

---

## Critical Components for Rust Port

### MUST IMPLEMENT (High Priority)

1. **Variant Struct** (equivalent to Java Variant.java)
   - All 30+ fields
   - Methods: `varType()`, `isGoodVar()`, `adjComplex()`, `isNoise()`

2. **BAM Parser** (using rust-htslib)
   - Query region from indexed BAM
   - Parse SAM records
   - Filter reads

3. **CIGAR Parser** ⭐
   - Parse CIGAR strings: M, I, D, N, S, H, etc.
   - Extract SNPs, indels, soft clips
   - Create Variation objects

4. **Variation Grouping**
   - HashMap<Position, Vec<Variation>>
   - Group by position and allele

5. **Statistics Calculator** ⭐
   - Position distribution (isAtLeastAt2Positions)
   - Quality distribution (hasAtLeast2DiffQualities)
   - Frequency calculation
   - Strand bias detection

6. **Soft-Clip Realigner** ⭐⭐
   - Extract soft-clipped sequences
   - Local sequence alignment (Smith-Waterman or similar)
   - Detect hidden indels

7. **Quality Filter**
   - Implement isGoodVar() logic
   - MSI detection
   - Strand bias calculation

8. **Output Formatter**
   - VCF-like text format
   - All fields from Variant struct

### NICE TO HAVE (Medium Priority)

1. Parallel region processing (rayon)
2. Multiple modes (Somatic, Amplicon, etc.)
3. Advanced filtering options
4. Performance optimizations

### NOT NEEDED FOR PHASE 1

1. Other modes (Somatic, Amplicon, Splicing)
2. Structural variant processing
3. Complex VCF handling
4. R integration for post-processing

---

## Complexity Assessment

### Simple (Easy to implement)
- Configuration parsing ✓
- BED file reading ✓
- Basic BAM reading ✓
- Output formatting ✓

### Medium (Straightforward)
- CIGAR parsing (well-defined spec)
- Read filtering
- Allele frequency calculation
- Strand bias detection

### Hard (Complex logic) ⭐⭐⭐
- **Soft-clip realignment** - Requires alignment algorithm, tricky edge cases
- **Statistics calculation** - Many distributions, many edge cases
- **Quality filtering** - Many interdependent thresholds
- **Complex variant handling** - Edge cases in coordinate representation

### Performance Critical ⭐
- Variation grouping (scale with depth)
- Statistics calculation (hot path)
- Output formatting (I/O bound)

---

## Files to Reference During Implementation

### Architecture
- `VARDICTJAVA_OVERVIEW.md` (this workspace) - High-level overview
- `SIMPLEMODE_PIPELINE.md` (this workspace) - Detailed pipeline walkthrough
- `JAVA_CODE_SNIPPETS.md` (this workspace) - Actual code examples

### Key Source Files to Study
1. **SimpleMode.java** (128 lines) - Entry point, pipeline creation
2. **Variant.java** (412 lines) - Data structure + key methods
3. **ToVarsBuilder.java** (1059 lines) - Statistics calculation LOGIC
4. **SimplePostProcessModule.java** (109 lines) - Quality filtering
5. **VariationRealigner.java** - Soft-clip realignment (key algorithm)

### Actual File Locations
```
VarDictJava/src/main/java/com/astrazeneca/vardict/
├── modes/SimpleMode.java
├── variations/Variant.java
├── modules/ToVarsBuilder.java
├── postprocessmodules/SimplePostProcessModule.java
└── modules/VariationRealigner.java
```

---

## Implementation Roadmap

### Phase 1: Core Simple Mode (Essential)
1. ✓ Configuration/CLI parser
2. ✓ Reference FASTA reader
3. ✓ BED region reader
4. BAM file parser (rust-htslib)
5. SAM record filtering
6. CIGAR parser
7. Variation tracking (HashMap)
8. Statistics calculator (ToVarsBuilder equivalent)
9. Soft-clip realigner (VariationRealigner equivalent)
10. Variant struct (Variant equivalent)
11. Quality filter (isGoodVar implementation)
12. Output formatter
13. Main pipeline orchestration

### Phase 2: Optimization
1. Parallel region processing
2. Memory-efficient data structures
3. SIMD for calculations
4. Cache optimization

### Phase 3: Features (Later)
1. Other modes (Somatic, Amplicon)
2. Additional filters
3. VCF output format
4. Advanced post-processing

---

## Performance Targets

Based on Java implementation (10x faster than Perl):
- Expected Rust performance: 15-20x vs Perl, 1.5-2x vs Java
- Bottleneck: I/O (BAM reading) - not CPU-bound
- Parallelizable: Region processing (good scaling)

### Optimization Opportunities
1. Reduce allocations in statistics calculation
2. Use SmallVec for position/quality vectors
3. String interning for alleles
4. Lazy quality distribution calculation
5. Cache-friendly data layout for variants

---

## Key Design Decisions for Rust

1. **Use rust-htslib** for BAM/SAM reading (battle-tested, performant)
2. **Struct-of-arrays** for Variant collections (better cache locality)
3. **Rayon** for parallel processing (simple API)
4. **HashMap** for variation grouping (not TreeMap, simpler)
5. **SmallVec** for small distributions (avoid heap)
6. **Arc<str>** for alleles (avoid duplication)

---

## Summary

VarDictJava is a **well-designed, modular variant caller** with clear separation of concerns:
- **Parsing layer** (SAM/CIGAR)
- **Processing layer** (realignment, statistics)
- **Output layer** (filtering, formatting)

The logic is **deterministic and testable**, with clear input/output for each module. The **hardest parts** are the statistics calculation and soft-clip realignment, but both are well-documented in the code.

For simple mode port, focus on:
1. Getting BAM/CIGAR parsing correct (unit test heavily)
2. Implementing exact statistics calculation (verify against Java)
3. Soft-clip realignment (most novel feature)
4. Quality filtering logic (many edge cases)

The rest is straightforward data transformation.

---

## Questions to Ask During Implementation

1. How to handle multiple variant alleles at same position?
2. How to calculate MSI accurately?
3. What alignment algorithm for soft-clip realignment?
4. How to handle complex variants (combine adjacent SNP+indel)?
5. How to determine optimal CIGAR when realigning?
6. What are the edge cases for isGoodVar()?

All answers are in the Java code - it's your source of truth.

---

**Analysis Date:** January 10, 2026  
**VarDictJava Version:** Final (Last maintained version)  
**Scope:** Simple Mode Single-Sample Variant Calling  
**Status:** Ready for Rust port implementation ✓

