# VarDict Java: Perl and R Script Analysis

## Overview
VarDict Java uses external Perl and R scripts for post-processing and statistical calculations. The original Perl-based VarDict is included as a Git submodule, and several critical scripts are called during the pipeline.

---

## External Scripts Used

### 1. **teststrandbias.R** (R Script - Strand Bias Testing)
**Location**: `VarDictJava/VarDict/teststrandbias.R`

**Purpose**: Performs strand bias statistical testing on variant calls

**Statistics Computed**:
- **Fisher Exact Test**: Tests strand bias by performing fisher.test on a 2x2 contingency table
  - Inputs: `d[i,10], d[i,11], d[i,12], d[i,13]` (forward/reverse counts for reference/variant)
  - Outputs: 
    - P-value (rounded to 5 decimal places)
    - Odds Ratio (rounded to 5 decimal places)

**Column Requirements**: 34, 36, or 38 columns
- 34: Standard bed files (pre-VarDictJava 1.5.5)
- 36: VarDictJava >= 1.5.5
- 38: Amplicon mode

**Java Equivalent**: Implemented in [src/main/java/com/astrazeneca/vardict/data/fishertest/FisherExact.java](VarDictJava/src/main/java/com/astrazeneca/vardict/data/fishertest/FisherExact.java)

---

### 2. **testsomatic.R** (R Script - Somatic Variant Testing)
**Location**: `VarDictJava/VarDict/testsomatic.R`

**Purpose**: Performs strand bias and somatic variant statistical testing on tumor-normal paired samples

**Statistics Computed**:
1. **Fisher Exact Test on Tumor Forward/Reverse**: 
   - Inputs: `d[i,10], d[i,11], d[i,12], d[i,13]`
   - Outputs: pvalues1, oddratio1

2. **Fisher Exact Test on Normal Forward/Reverse**:
   - Inputs: `d[i,28], d[i,29], d[i,30], d[i,31]`
   - Outputs: pvalues2, oddratio2

3. **Fisher Exact Test on Somatic Status** (one-sided):
   - Inputs: Variant coverage, reference coverage for tumor and normal
   - Formula: `tref = max(0, d[i,8] - d[i,9])` and `rref = max(0, d[i,26] - d[i,27])`
   - Tests both "greater" and "less" alternatives, uses the one with lower p-value
   - Outputs: pvalues, oddratio

**Column Requirements**: >= 48 columns

**Java Equivalent**: Implemented in [src/main/java/com/astrazeneca/vardict/data/fishertest/FisherExact.java](VarDictJava/src/main/java/com/astrazeneca/vardict/data/fishertest/FisherExact.java)

---

### 3. **var2vcf_valid.pl** (Perl Script - VCF Conversion for Single Sample)
**Location**: `VarDictJava/VarDict/var2vcf_valid.pl`

**Purpose**: Converts VarDict output to VCF format for single sample/germline variants

**Functions**:
- Filters variants based on quality criteria
- Converts variant output to VCF v4.2 format
- Applies filtering thresholds (quality, depth, frequency, etc.)

**Configuration Parameters**:
- `$TotalDepth`: Minimum total depth (default: 3)
- `$VarDepth`: Minimum variant depth (default: 2)
- `$Freq`: Minimum frequency (default: 0.02)
- `$Pmean`: Mean position (default: 8)
- `$qmean`: Base quality mean (default: 22.5)
- `$Qmean`: Mapping quality mean (default: 10)
- `$GTFreq`: Genotype frequency (default: 0.2)
- `$SN`: Signal to Noise ratio (default: 1.5)

**Java Status**: **NOT PORTED** - Currently requires Perl

---

### 4. **var2vcf_paired.pl** (Perl Script - VCF Conversion for Tumor-Normal Pairs)
**Location**: `VarDictJava/VarDict/var2vcf_paired.pl`

**Purpose**: Converts VarDict output to VCF format for tumor-normal paired samples

**Functions**:
- Filters somatic and germline variants
- Applies variant classification logic (somatic, germline, LOH, etc.)
- Converts variant output to VCF v4.2 format

**Java Status**: **NOT PORTED** - Currently requires Perl

---

## Statistical Calculations Breakdown

### Fisher Exact Test (Primary Statistical Test)
**Implementation**: [FisherExact.java](VarDictJava/src/main/java/com/astrazeneca/vardict/data/fishertest/FisherExact.java)

**What it calculates**:
- Non-central hypergeometric distribution parameters
- Maximum Likelihood Estimation (MLE) of odds ratio
- P-values for strand bias testing
- Support: one-sided (less/greater) and two-sided tests

**Key Methods**:
- `FisherExact(int refFwd, int refRev, int altFwd, int altRev)`: Constructor
- `getPValue()`: Returns two-sided p-value
- `getPValueLess()`: Returns one-sided p-value (less)
- `getPValueGreater()`: Returns one-sided p-value (greater)
- `getOddRatio()`: Returns odds ratio as string
- Uses Apache Commons Math3's `HypergeometricDistribution`

**Replacement Status**: ✅ ALREADY PORTED TO JAVA - Can be replicated in Rust

---

### Basic Statistical Calculations (Perl Module)
**Location**: `VarDictJava/VarDict/Stat/Basic.pm`

**Functions Available** (for reference):
- `mean()` - Calculate mean
- `sum()` - Calculate sum
- `min()` - Find minimum
- `max()` - Find maximum
- `median()` - Calculate median
- `var()` - Calculate variance
- `std()` - Calculate standard deviation
- `mad()` - Mean/Median Absolute Deviation
- `rstd()` - Robust standard deviation (MAD * 1.4826)
- `prctile()` - Calculate percentile
- `iqr()` - Calculate interquartile range
- `zscore()` - Calculate z-scores
- `standardize()` - Standardize values

**Java Status**: **BUILT-IN OR AVAILABLE** - These are standard statistics available in Java/Rust

---

## Java Fisher Exact Test Usage

The Java implementation is used in three output variants:

1. **SimpleOutputVariant.java** - Single sample output
   ```java
   FisherExact fisher = new FisherExact(
       variant.refForwardCoverage, 
       variant.refReverseCoverage,
       variant.variantForwardCoverage, 
       variant.variantReverseCoverage
   );
   this.pvalue = fisher.getPValue();
   this.oddratio = fisher.getOddRatio();
   ```

2. **SomaticOutputVariant.java** - Tumor-normal paired samples
   ```java
   // Tumor strand bias test
   FisherExact fisher = new FisherExact(
       tumorVariant.refForwardCoverage, 
       tumorVariant.refReverseCoverage,
       tumorVariant.variantForwardCoverage, 
       tumorVariant.variantReverseCoverage
   );
   
   // Normal strand bias test
   fisher = new FisherExact(
       normalVariant.refForwardCoverage, 
       normalVariant.refReverseCoverage,
       normalVariant.variantForwardCoverage, 
       normalVariant.variantReverseCoverage
   );
   
   // Somatic comparison test
   fisher = new FisherExact(
       this.var1variantCoverage, tref, 
       this.var2variantCoverage, rref
   );
   ```

3. **AmpliconOutputVariant.java** - Amplicon mode output

---

## Configuration Flag: --fisher

When the `--fisher` flag is enabled in VarDictJava:
- The Java implementation of Fisher Exact Test is used **INSTEAD OF** calling R scripts
- Files: `teststrandbias.R` and `testsomatic.R` are **NOT NEEDED**
- This is an "Experimental feature" per the code comments

**CmdParser Reference**:
```java
config.fisher = cmd.hasOption("fisher");
```

---

## Summary for Rust Implementation

### ✅ Already Ported (Can replicate directly):
1. **Fisher Exact Test** - Complete Java implementation exists
   - Uses non-central hypergeometric distribution
   - Calculates odds ratio via MLE
   - Returns p-values (one-sided and two-sided)
   - **Action**: Implement in Rust using statistical libraries

### ⚠️ Partially Ported (Need Perl):
1. **var2vcf_valid.pl** - VCF output formatting
2. **var2vcf_paired.pl** - VCF output formatting for paired samples

### ✅ Not Needed for Core Variant Detection:
1. **teststrandbias.R** - Replaced by Java FisherExact when `--fisher` flag is used
2. **testsomatic.R** - Replaced by Java FisherExact when `--fisher` flag is used
3. **Stat/Basic.pm** - Basic math functions available in all languages

---

## Recommendations for Rust Port

### Priority 1 - Required for Core Functionality:
- [ ] Implement Fisher Exact Test with hypergeometric distribution
- [ ] Use a Rust statistical library (e.g., `statrs`, `ndarray-stats`, or similar)

### Priority 2 - Post-Processing:
- [ ] Implement VCF output formatting (var2vcf_valid.pl logic)
- [ ] Implement VCF output for paired samples (var2vcf_paired.pl logic)
- [ ] Consider optional integration with external R if needed

### Priority 3 - Optional:
- [ ] Basic statistics module (mean, std, median, etc.) - Can use existing Rust crates
- [ ] Custom Perl script support if backward compatibility is needed
