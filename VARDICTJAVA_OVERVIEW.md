# VarDict Java - Architecture & Logic Overview

## Summary
VarDict is a sensitive variant discovery program from next-generation sequencing (NGS) data. It's written in Java (10x faster than the original Perl version) and implements features like amplicon bias awareness, soft-clipped read realignment, and good scalability.

**Current Status:** Final Version - No longer maintained (Java port)

## Project Structure

### Directory Layout
```
VarDictJava/
├── src/main/java/com/astrazeneca/vardict/
│   ├── Main.java                 # Entry point
│   ├── Configuration.java        # Config management
│   ├── CmdParser.java           # Command-line parser
│   ├── VarDictLauncher.java     # Launcher
│   ├── Utils.java               # Utility functions
│   ├── RegionBuilder.java       # BED region parser
│   ├── modes/                   # Different analysis modes
│   │   ├── SimpleMode.java      # ⭐ FOCUS: Simple single-sample mode
│   │   ├── AbstractMode.java    # Base mode class
│   │   ├── AmpliconMode.java    # Amplicon mode
│   │   ├── SomaticMode.java     # Paired sample mode
│   │   └── SplicingMode.java    # Splicing mode
│   ├── data/                    # Data structures
│   │   ├── Reference.java       # Reference genome
│   │   ├── Region.java          # BED region info
│   │   ├── ReferenceResource.java # Resource management
│   │   ├── patterns.rs          # Regex patterns
│   │   └── scopedata/           # Pipeline scope data
│   │       ├── Scope.java       # Common pipeline scope
│   │       ├── InitialData.java # Initial BAM data
│   │       ├── AlignedVarsData.java # Aligned variants
│   │       └── GlobalReadOnlyScope.java # Global settings
│   ├── modules/                 # Pipeline modules
│   │   ├── SAMFileParser.java   # BAM file reader
│   │   ├── RecordPreprocessor.java # SAM record preprocessing
│   │   ├── CigarParser.java     # CIGAR string parser
│   │   ├── CigarModifier.java   # CIGAR modification
│   │   ├── VariationRealigner.java # Realign soft clips
│   │   ├── ToVarsBuilder.java   # Build variant objects ⭐
│   │   └── StructuralVariantsProcessor.java # SV processing
│   ├── variations/              # Variant data structures
│   │   ├── Variant.java         # ⭐ Main variant class
│   │   ├── Vars.java            # Collection at a position
│   │   ├── Variation.java       # Low-level variation
│   │   ├── VariationUtils.java  # Variation utilities
│   │   └── others...
│   ├── postprocessmodules/      # Post-processing
│   │   ├── SimplePostProcessModule.java # ⭐ Output variants
│   │   ├── SomaticPostProcessModule.java
│   │   └── AmpliconPostProcessModule.java
│   ├── printers/                # Output formatters
│   │   ├── VariantPrinter.java
│   │   └── SimpleOutputVariant.java
│   └── collection/              # Custom collections
```

## Pipeline Architecture (Simple Mode)

### Overall Data Flow
```
BED File + BAM File + Reference FASTA
         ↓
   Configuration
         ↓
   VarDictLauncher
         ↓
   SimpleMode.notParallel() or .parallel()
         ↓
  ┌─────────────────────────────────────┐
  │  For each Region (Sequential/Parallel)
  └─────────────────────────────────────┘
         ↓
  ┌─────────────────────────────────────────────────────────┐
  │           processBamInPipeline(region)                   │
  │                                                           │
  │  1. SAMFileParser → Read BAM records for region          │
  │  2. RecordPreprocessor → Filter/validate reads          │
  │  3. CigarParser → Parse CIGAR strings                   │
  │  4. CigarModifier → Adjust CIGAR for soft clips         │
  │  5. VariationRealigner → Realign clipped bases          │
  │  6. ToVarsBuilder → Create Variant objects ⭐           │
  │  7. SimplePostProcessModule → Filter variants           │
  │  8. VariantPrinter → Output VCF/text format             │
  └─────────────────────────────────────────────────────────┘
         ↓
   VCF/Text Output
```

## Key Classes Explained

### 1. **SimpleMode.java** - Entry Point for Simple Single-Sample Analysis
**Purpose:** Orchestrates the variant calling pipeline for single samples

**Key Methods:**
- `notParallel()` - Processes regions sequentially
- `parallel()` - Creates worker threads for parallel region processing
- `processBamInPipeline(region, printer)` - Main pipeline execution

**Flow:**
1. For each region in segments:
2. Create a `Scope<InitialData>` with BAM file, region info, reference
3. Chain modules: Parser → Realigner → ToVarsBuilder → PostProcessor
4. Print results

```java
// Pseudocode pipeline
Scope<InitialData> → SAMFileParser 
                   → RecordPreprocessor 
                   → CigarParser 
                   → VariationRealigner 
                   → ToVarsBuilder 
                   → SimplePostProcessModule 
                   → Output
```

### 2. **Scope<T>** - Pipeline Scope Container
**Purpose:** Carries data through the pipeline stages

**Contains:**
- `bam` - BAM file path
- `region` - Current BED region
- `regionRef` - Reference sequence for region
- `referenceResource` - Reference genome
- `maxReadLength` - Longest read length
- `splice` - Splice sites
- `out` - Output printer
- `data` - Stage-specific data (generic T)

**Generic Progression:**
```
Scope<InitialData> 
  → Scope<RealignedVariationData> 
  → Scope<AlignedVarsData>
```

### 3. **Variant.java** - Variant Data Structure
**Purpose:** Holds all variant information for output

**Key Fields:**
- `descriptionString` - Variant description (SNP letter, +seq for ins, -count for del, etc.)
- `positionCoverage` - Depth at variant position
- `varsCountOnForward/Reverse` - Strand counts
- `frequency` - Allele frequency
- `meanPosition` - Mean position in read
- `meanQuality` - Mean base quality
- `meanMappingQuality` - Mean mapping quality
- `strandBiasFlag` - Strand bias (0, 1, or 2)
- `startPosition` / `endPosition` - Variant coordinates
- `refallele` / `varallele` - Reference and variant alleles
- `vartype` - Variant type (SNP, Insertion, Deletion, Complex)
- `isGoodVar()` - Quality filter check

**Description Format:**
1. Single letter (SNP): "A", "T", "G", "C"
2. Insertion: "+ACGT"
3. Deletion: "-5" (5 bases deleted)
4. Complex indel: "ACGT#-3" (insertion followed by deletion)
5. Followed by insertion: "^ACGT"
6. Followed by deletion: "^5"

### 4. **Vars.java** - Variants at a Position
**Purpose:** Collection of variants at a single genomic position

**Contains:**
- `variants` - List of `Variant` objects at position
- `referenceVariant` - Reference "variant" (for pileup mode)
- `sv` - Structural variant flags

### 5. **ToVarsBuilder.java** - Variant Creation ⭐⭐⭐
**Purpose:** Converts low-level variations into high-quality variant objects

**Key Operations:**
1. Collects all variations from realignment step
2. Groups by position
3. Calculates statistics:
   - Allele frequencies
   - Mean positions in reads
   - Base quality stats
   - Mapping quality stats
   - Strand bias
4. Filters variants based on:
   - Frequency threshold
   - Quality thresholds
   - Strand bias
   - MSI (microsatellite) detection
5. Creates final `Variant` objects

**Input:** `Scope<RealignedVariationData>`
**Output:** `Scope<AlignedVarsData>` (map of position → Vars)

### 6. **SimplePostProcessModule.java** - Output Filtering ⭐
**Purpose:** Final filtering and output formatting for simple mode

**Process:**
1. Iterate through aligned variants map
2. For each position:
   - Skip if outside region (unless SVs present)
   - Filter reference calls (unless pileup mode)
   - Check variant quality: `isGoodVar()`
   - Handle complex variants: `adjComplex()`
   - Apply CRISPR cut site filtering (if enabled)
3. Create `SimpleOutputVariant` objects
4. Send to printer

**Key Filters:**
- Variant must be in region bounds
- Reference allele can't contain 'N'
- Must pass `isGoodVar()` quality check
- Complex variants adjusted for CRISPR sites

### 7. **SAMFileParser.java** - BAM Input
**Purpose:** Reads BAM file for a region

**Steps:**
1. Open BAM file
2. Query region using SAM/BAM index
3. Parse SAM records
4. Convert to internal representation
5. Pass to next module

### 8. **VariationRealigner.java** - Soft Clip Realignment
**Purpose:** Rescue long indels hidden in soft-clipped sequences

**Key Feature:** One of VarDict's novel features - realigns soft-clipped reads to find hidden variants

**Process:**
1. Extract soft-clipped sequences
2. Realign to reference
3. Find indels that would explain the clipping
4. Add as new variations

### 9. **CigarParser & CigarModifier** - CIGAR Handling
**Purpose:** Parse and modify CIGAR strings for variant detection

**CIGAR Operations:**
- M = alignment match (can be SNP, match, or mismatch)
- I = insertion to reference
- D = deletion from reference
- N = skipped region
- S = soft clipping
- H = hard clipping

## Configuration (Configuration.java)

**Key Settings:**
- `bam` - BAM file path(s)
- `regions` - BED file regions
- `ref` - Reference FASTA file
- `minFreq` - Minimum allele frequency (default: 0.02)
- `minQuality` - Minimum base quality
- `minMappingQuality` - Minimum MAPQ
- `minAltCount` - Minimum variant count
- `doPileup` - Output reference calls (default: false)
- `crisprCuttingSite` - CRISPR cut site adjustment
- `threads` - Number of parallel threads
- `printerType` - Output format (VCF, text, etc.)

## Important Features for Simple Mode

### 1. **Soft-Clipped Read Realignment**
- Rescue variants hidden in soft clips
- Better indel detection

### 2. **Strand Bias Detection**
- Flag variants with strong strand bias
- Values: 0 (no bias), 1 (weak), 2 (strong)

### 3. **Mean Position & Quality Variance**
- Track variant position across reads
- Track quality variation (multimodal)
- Flags artifacts

### 4. **Microsatellite Detection**
- Identify MSI regions (repeating sequences)
- Adjust frequency thresholds in MSI areas

### 5. **Complex Variant Handling**
- Combine adjacent SNPs + indels
- Adjust coordinates for better representation

### 6. **CRISPR Cut Site Awareness**
- Mark variants near CRISPR cut sites
- Useful for genome-edited sample analysis

## Output Format

Default is **VCF-like text format** with fields:
```
CHROM POS REF ALT DEPTH ALT_COUNT ALT_FREQ REFALLELE VARALLELE VARTYPE
      [additional fields for strand bias, quality, etc.]
```

## Execution Modes

### Sequential (not parallel)
```
Thread 1: Region 1 → Region 2 → Region 3 ...
```
**Pros:** Simple, low memory
**Cons:** Slow on multi-core systems

### Parallel
```
Thread 1: Region 1
Thread 2: Region 2
Thread 3: Region 3
... (executor pool)
Output: Queue-based, maintains input order
```
**Pros:** Fast on multi-core
**Cons:** Higher memory usage

## Key Algorithms

### 1. Variant Calling
- Count base occurrences at each position
- Group by identity (SNP) or offset (indel)
- Calculate statistics from read attributes
- Filter by frequency and quality

### 2. Indel Detection
- Parse CIGAR string for I/D operations
- Extract indel sequences
- Realign soft clips for hidden indels
- Combine adjacent indels (complex)

### 3. Quality Filtering
- Position in read (should be middle, not ends)
- Base quality (typically ≥ Q20)
- Mapping quality (typically ≥ Q20)
- Strand bias (low threshold in simple mode)
- Variant frequency (typically ≥ 2%)

### 4. Realignment
- Local sequence alignment
- Find optimal indel representation
- 3' shifting for deletions (right-normalization)

## Performance Characteristics

- **Java vs Perl:** ~10x faster
- **I/O Bound:** BAM file reading is bottleneck
- **CPU:** Region processing can be parallelized
- **Memory:** ~1-2GB for typical operations
- **Scalability:** Good for whole genome (multiple workers)

## Next Steps for Rust Port

### Phase 1: Core Infrastructure (Simple Mode only)
1. ✅ Reference FASTA parsing
2. ✅ BED region reading
3. ✅ Configuration/CLI parsing
4. ✅ BAM file reading (using rust-htslib)
5. SAM record parsing & filtering
6. CIGAR parsing
7. Variation tracking (hash maps)
8. Soft clip realignment
9. Variant creation (statistics calculation)
10. Post-processing & filtering
11. VCF output

### Phase 2: Optimization
- Parallel region processing
- Memory-efficient data structures
- SIMD for statistics calculation
- Custom allocator for frequency maps

### Phase 3: Features (Later)
- Other modes (Somatic, Amplicon, etc.)
- Advanced filtering
- Additional output formats

## Key Files for Understanding

1. **src/main/java/com/astrazeneca/vardict/modes/SimpleMode.java** - Pipeline orchestration
2. **src/main/java/com/astrazeneca/vardict/variations/Variant.java** - Variant data structure
3. **src/main/java/com/astrazeneca/vardict/modules/ToVarsBuilder.java** - Variant creation logic
4. **src/main/java/com/astrazeneca/vardict/postprocessmodules/SimplePostProcessModule.java** - Output filtering
5. **src/main/java/com/astrazeneca/vardict/modules/VariationRealigner.java** - Key algorithm
6. **src/main/java/com/astrazeneca/vardict/data/scopedata/Scope.java** - Pipeline data flow

---

**Last Updated:** January 10, 2026
**Based on:** VarDictJava source code exploration
**Focus:** Simple Mode Single-Sample Variant Calling
