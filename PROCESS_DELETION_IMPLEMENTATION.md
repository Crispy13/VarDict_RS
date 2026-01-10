# Implementation: process_deletion() in Rust

## Summary
Successfully implemented the `process_deletion()` function in Rust for VarDict's CIGAR parser, using the `VarDesc` enum instead of string descriptions.

## File Modified
`/home/eck/workspace/vardict_rs/src/mods/cigar_parser.rs` (lines 622-850)

## Key Implementation Details

### VarDesc Enum Usage
Instead of Java's string-based deletion description (`"-5"` for 5-base deletion), the Rust implementation uses:
```rust
pub(crate) enum VarDesc {
    Del { len: u32 }  // Type-safe deletion representation
}
```

### Algorithm Implementation

The function handles three main deletion scenarios:

#### 1. Complex Indel Pattern (D + M + I/D)
- Scans through a matched segment after deletion
- Detects mismatches and tracks position variance
- Updates both read and reference position offsets

#### 2. Deletion Followed by Insertion (D + I)
- Appends insertion sequence with `^` marker
- Tracks quality of inserted segment
- Skips the next insertion in CIGAR processing

#### 3. Deletion Followed by Match (D + M)
- Scans through matched bases after deletion
- Looks for good quality bases to anchor the deletion
- Appends matching sequence if found

### Key Features Implemented

1. **Intron-Adjacent Deletion Detection**
   - Skips deletions next to introns (RNA-seq artifacts)
   - Early return if detected

2. **Quality-Based Variant Tracking**
   ```rust
   var.mean_pos += tp as f64;        // Position in read
   var.mean_qual += tmpq;             // Quality average
   var.mean_mapq += mapq as f64;      // Mapping quality
   var.nm += (nm - nmoff) as f64;     // Mismatch count
   ```

3. **Position Variance Calculation**
   - Calculates position in read relative to ends
   - `tp = min(readPos, totalLen - readPos)`
   - Used for quality filtering later

4. **Reference Coverage Tracking**
   - Increments coverage for all deleted bases
   - Important for allele frequency calculations
   - Loop: `for i in 0..cigarElementLength`

5. **Complex Indel Support**
   - Handles adjacent insertions/deletions
   - Tracks mismatch offsets separately
   - Skips processed segments correctly

### Parameter Explanpts

| Parameter | Purpose |
|-----------|---------|
| `query_sequence` | Read bases |
| `mapq` | Mapping quality of the read |
| `contig_ref_seq` | Reference sequence for region |
| `query_quality` | Base quality scores (ASCII - 33) |
| `nm` | Number of mismatches in read |
| `is_reverse` | Strand direction (true = reverse) |
| `read_len_including_match_ins` | Total read length |
| `ci` | Current CIGAR element index |

### Configuration Parameters Used

- `conf.vext` - Window size for looking ahead (default: ~10bp)
- `conf.goodq` - Minimum quality threshold (default: 20)
- `conf.perform_local_realignment` - Enable local realignment (soft-clip detection)
- `conf.goodq` - Quality threshold for high-quality reads

### Data Structures Modified

1. **Variant struct fields updated:**
   - `alt_depth` - Count incremented
   - `alt_depth_fwd/alt_depth_rev` - Strand-specific count
   - `mean_pos` - Position distribution
   - `mean_qual` - Quality distribution
   - `mean_mapq` - Mapping quality
   - `nm` - Mismatch count
   - `high_qual_read_cnt/low_qual_read_cnt` - Quality-based counts

2. **CigarParser fields updated:**
   - `non_insertion_vars` - Map of deletion variants
   - `ref_coverage` - Coverage for deleted positions
   - `read_pos_including_softclip` - Read position including soft clips
   - `read_pos_excluding_softclip` - Read position excluding soft clips
   - `start` - Reference position
   - `offset` - Offset for complex indels

## Algorithm Comparison with Java

### Java Approach
```java
// String-based deletion description
StringBuilder descStringOfDeletedElement = new StringBuilder("-").append(cigarElementLength);
// Appended to nonInsertionVariants map
Variation hv = getVariation(nonInsertionVariants, start, descStringOfDeletedElement.toString());
```

### Rust Approach
```rust
// Type-safe enum-based deletion description
let del_desc = VarDesc::Del { len: self.cigar_len };
// Used as HashMap key
let var: &mut Variant = get_variants_from_map(&mut self.non_insertion_vars, self.start, &del_desc);
```

**Advantages:**
- Type-safe (no string parsing errors)
- Efficient hashing (u32 vs String)
- Compile-time correctness
- Better performance (no allocations)

## Edge Cases Handled

1. **Out-of-bounds quality access**
   - Uses `saturating_sub()` for safe indexing
   - Checks bounds before array access

2. **Intron-adjacent indels**
   - Early exit for RNA-seq artifacts
   - Uses `skip_indel_next_to_intron()` helper

3. **Complex quality calculations**
   - Empty quality vector handled (prevents division by zero)
   - Best-of-two-qualities for anchoring bases

4. **Reference position adjustments**
   - Accounts for multiple indel patterns
   - Updates both read and reference offsets correctly

5. **Coverage tracking**
   - Increments for all deleted bases in reference
   - Important for copy number variation detection

## Testing Recommendations

1. **Unit tests needed:**
   - Simple 5-base deletion
   - Complex deletion + match pattern
   - Deletion + insertion pattern
   - Deletion adjacent to intron (skip)
   - Low-quality bases after deletion

2. **Integration tests:**
   - Compare output with Java VarDict
   - Verify coverage calculations
   - Check position variance calculations
   - Validate quality distributions

3. **Comparison points:**
   ```
   Java VarDict | Rust VarDict
   "-5" string  | VarDesc::Del { len: 5 }
   Variant.pp   | (position variance tracking - TODO)
   Variant.pq   | (quality variance tracking - TODO)
   ```

## TODO / Future Improvements

1. **Position Variance Tracking (pstd flag)**
   - Add to Variant struct: `pp: Option<usize>` and `pstd: bool`
   - Track previous position to detect 2+ different positions

2. **Quality Variance Tracking (qstd flag)**
   - Add to Variant struct: `pq: Option<f64>` and `qstd: bool`
   - Track previous quality to detect 2+ different qualities

3. **Performance Optimization**
   - Use SmallVec for quality_of_segment (avoid heap allocation)
   - Consider SIMD for quality averaging

4. **Extended Deletion Handling**
   - Support for 3-base or 4-base matching patterns
   - Better handling of homopolymer runs

## Verification Status

✅ **Compilation:** Passes (warnings only for pre-existing issues)  
✅ **Type Safety:** Full enum usage for VarDesc  
✅ **Error Handling:** Proper Result returns and error propagation  
⏳ **Testing:** Needs unit and integration tests  
⏳ **Variance Tracking:** Needs pstd/qstd implementation  

## Files Referenced

**Java Implementation:**
- `VarDictJava/src/main/java/com/astrazeneca/vardict/modules/CigarParser.java` (lines 671-900)

**Rust Implementation:**
- `src/mods/cigar_parser.rs` (lines 622-850)
- `src/variants/variants.rs` (VarDesc enum)
- `src/variants/var_utils.rs` (get_variants_from_map function)

**Documentation Used:**
- SIMPLEMODE_PIPELINE.md (STEP 3: CIGAR parsing)
- JAVA_CODE_SNIPPETS.md (CigarParser logic)

---

**Implementation Date:** January 10, 2026  
**Status:** Complete (ready for testing)  
**Next Step:** Implement process_insertion() and unit tests

