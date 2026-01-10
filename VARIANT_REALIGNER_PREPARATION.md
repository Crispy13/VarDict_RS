# VariantRealigner Preparation

## Overview
- **Java file**: `VarDictJava/src/main/java/com/astrazeneca/vardict/modules/VariationRealigner.java` (2744 lines)
- **Rust file**: `src/mods/variant_realigner.rs` (110 lines, mostly stubs)
- **Purpose**: Realign variations (softclips, indels, long insertions) after initial variant calling

## Key Responsibilities

### Main Process Flow (process method)
1. Filter all SV structures (filterAllSVStructures)
2. Adjust MNP (adjustMNP)
3. Realign indels (realignIndels) - if perform_local_realignment enabled

### Major Methods Needed
1. **filterAllSVStructures()** - Filter possible structural variants (8 types: inv3, rinv3, inv5, rinv5, del, rdel, dup, rdup)
2. **filterSV(List<Sclip>)** - Filter individual SV list by checking clusters
3. **adjustMNP()** - Adjust multi-nucleotide polymorphisms
4. **realignIndels()** - Main realignment logic for indels
5. **realign_del()** - Realign deletions (partial stub exists)
6. **realign_ins()** - Realign insertions

### Helper Methods (with tests in Java)
1. **islowcomplexseq(String)** - Static method to detect low-complexity sequences
   - Test cases: "AAAAAAAAA" (true), "ATATATATATAT" (true), "CCCCCCCCGA" (true)
   - "ACGTACGTACGT" (false), "CCGTAACGGGGT" (false)

2. **ismatch(seq1, seq2, threshold)** - Instance method to check if sequences match within threshold
   - Test cases provided in VariationRealignerTest.java

3. **find35match(seq5, seq3)** - Find 3'/5' end matches between two sequences
   - Returns Match35 object with (matched5end, matched3End, maxMatchedLength)
   - Test data provided with expected results

## Data Structures Used
- `svStructures`: Contains all SV types (inv3, rinv3, inv5, rinv5, del, rdel, dup, rdup, fus, rfus)
- `Sclip`: Soft-clip data structure
- `Cluster`: Result of clustering analysis
- `Match35`: Result of 3'/5' match finding
- `nonInsertionVariants`: Map of position -> variant descriptions
- `insertionVariants`: Map of position -> variant descriptions
- `refCoverage`: Map of position -> coverage count
- `mnp`: Multi-nucleotide polymorphism data

## Java Unit Tests Available
✅ **testIsLowComplexSeq()** - 5 test cases
✅ **testIsMatch()** - 7 test cases  
✅ **testFind35Match()** - 5 data provider test cases

## Next Steps (When Starting Implementation)
1. **First**: Port utility helper functions (islowcomplexseq, ismatch, find35match)
2. **Second**: Port filterAllSVStructures and related SV filtering
3. **Third**: Port adjustMNP
4. **Fourth**: Port realignIndels and sub-methods

## Complexity Assessment
- **Large file**: 2744 lines with many nested structures
- **Many dependencies**: Uses Cluster, Sclip, SVStructures, VariationMap
- **Unit tests available**: Good coverage for helper methods
- **Recommendation**: Port in logical chunks, starting with tested utility functions
