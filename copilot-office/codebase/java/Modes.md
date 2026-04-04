# Modes

**Source**: `VarDictJava/src/main/java/com/astrazeneca/vardict/modes/`  
**Files**: `AbstractMode.java` (145 LOC), `SimpleMode.java` (118 LOC), `SomaticMode.java` (148 LOC), `AmpliconMode.java` (146 LOC), `SplicingMode.java` (112 LOC)  
**Total LOC**: ~669  
**Rust counterpart**: `src/mods/` (pipeline orchestration)  
**Status**: complete

---

## Overview

The Modes module is the **top-level orchestrator** of the VarDict pipeline. It determines which variant calling strategy to execute (single-sample, tumor/normal, amplicon, or splicing) and manages the lifecycle of all pipeline stages.

### Mode Selection Logic (VarDictLauncher.java:52-66)

```
if outputSplicing → SplicingMode
else if regionOfInterest != null OR ampliconBasedCalling == null:
    if hasBam2 → SomaticMode
    else → SimpleMode
else → AmpliconMode
```

After construction, the mode is stored in `GlobalReadOnlyScope.mode` via `setMode()`. This singleton is called back by `StructuralVariantsProcessor` and `VariationRealigner` and `SomaticPostProcessModule` via `getMode().partialPipeline()` or `getMode().pipeline()` for re-entrant sub-pipelines.

### Thread Model

| Scenario | Executor | Behavior |
|----------|----------|----------|
| `threads == 1` | `DirectThreadExecutor` | All work on calling thread, synchronous |
| `threads > 1` | `Executors.newFixedThreadPool(threads)` | Regions submitted to pool; results collected via `LinkedBlockingQueue(capacity=10)` with sentinel `LAST_SIGNAL_FUTURE` |
| Re-entrant calls (partialPipeline, pipeline from PostProcess) | Always `DirectThreadExecutor` | Always synchronous on calling thread |

### Key Data Structures

- **`Scope<T>`**: Immutable data carrier passed between pipeline stages. Contains `bam` (String path), `region`, `regionRef`, `referenceResource`, `maxReadLength`, `splice` (Set<String>), `out` (VariantPrinter), and `data` (the generic payload for current stage).
- **`DirectThreadExecutor`**: `Executor` that runs `command.run()` on the calling thread — makes `CompletableFuture.supplyAsync()` effectively synchronous.
- **`LAST_SIGNAL_FUTURE`**: `CompletableFuture.completedFuture(null)` used as a sentinel to signal the consumer thread that all tasks are done.

---

## Method Inventory

### AbstractMode.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `AbstractMode(List<List<Region>>, ReferenceResource)` | L30-L33 | yes | Constructor: stores segments and referenceResource |
| `pipeline(Scope<InitialData>, Executor)` | L42-L57 | yes | Full pipeline: SAMFileParser → CigarParser(false) → VariationRealigner → StructuralVariantsProcessor → ToVarsBuilder |
| `stopVardictWithException(Region, Throwable)` | L64-L69 | yes | Prints error with region coordinates, calls `System.exit(1)` |
| `partialPipeline(Scope<InitialData>, Executor)` | L79-L83 | yes | Partial pipeline: SAMFileParser → CigarParser(true) |
| `splicingPipeline(Scope<InitialData>, Executor)` | L93-L97 | yes | Splicing pipeline: SAMFileParser → CigarParser(false) |
| `notParallel()` | L99 | yes | Abstract: sequential mode |
| `parallel()` | L101-L103 | yes | Calls `createParallelMode().process()` |
| `createParallelMode()` | L105 | yes | Abstract: creates parallel mode worker |
| `AbstractParallelMode.process()` | L114-L134 | yes | Parallel execution loop: consumer thread takes from queue, prints results |
| `AbstractParallelMode.produceTasks()` | L136 | yes | Abstract: producer thread submits tasks to queue |
| `printHeader()` | L141 | yes | Abstract: prints TSV header |
| `tryToGetReference(Region)` | L143-L150 | yes | Gets reference sequence, catches exception with `stopVardictWithException` |

### SimpleMode.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `SimpleMode(List<List<Region>>, ReferenceResource)` | L32-L35 | yes | Constructor; calls `printHeader()` |
| `notParallel()` | L41-L47 | yes | Iterates all segments/regions, calls `processBamInPipeline` |
| `createParallelMode()` | L53-L63 | yes | Creates anonymous `AbstractParallelMode` that submits `VardictWorker` per region |
| `VardictWorker(Region)` | L69-L83 | yes | `Callable<OutputStream>`: creates `ByteArrayOutputStream`, runs `processBamInPipeline`, returns stream |
| `processBamInPipeline(Region, VariantPrinter)` | L90-L104 | yes | Creates `Scope<InitialData>`, runs `pipeline()`, chains `SimplePostProcessModule`, joins |
| `printHeader()` | L106-L118 | yes | Prints 36-column header (+ optional CRISPR column) |

### SomaticMode.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `SomaticMode(List<List<Region>>, ReferenceResource)` | L33-L36 | yes | Constructor; calls `printHeader()` |
| `notParallel()` | L42-L51 | yes | Iterates segments/regions; creates splice set and ref per region; calls `processBothBamsInPipeline` |
| `createParallelMode()` | L57-L70 | yes | Creates anonymous `AbstractParallelMode` that submits `SomdictWorker` per region |
| `SomdictWorker(Region, Set<String>, Reference)` | L76-L95 | yes | `Callable<OutputStream>`: creates stream, runs `processBothBamsInPipeline`, returns stream |
| `processBothBamsInPipeline(VariantPrinter, Region, Set<String>, Reference)` | L97-L125 | yes | Runs pipeline on bam1, joins, runs pipeline on bam2 with bam1's maxReadLength, chains `SomaticPostProcessModule.thenAcceptBoth`, joins |
| `printHeader()` | L127-L140 | yes | Prints somatic header (duplicated tumor/normal column sets) |

### AmpliconMode.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `AmpliconMode(List<List<Region>>, ReferenceResource)` | L36-L39 | yes | Constructor; calls `printHeader()` |
| `notParallel()` | L45-L78 | yes | Iterates segments; builds position→amplicon map; runs pipeline per amplicon region; collects vars; calls `AmpliconPostProcessModule` |
| `createParallelMode()` | L84-L131 | yes | Creates anonymous parallel mode: builds pos map, submits pipelines per region to executor, collects results, submits AmpliconPostProcessModule as async task |
| `printHeader()` | L133-L142 | yes | Prints amplicon header (includes GoodVarCount, TotalVarCount, Nocov, Ampflag) |

### SplicingMode.java

| Method | Lines | Analyzed? | Summary |
|--------|-------|-----------|---------|
| `SplicingMode(List<List<Region>>, ReferenceResource)` | L29-L32 | yes | Constructor; calls `printHeader()` |
| `notParallel()` | L38-L44 | yes | Iterates segments/regions, calls `processRegion` |
| `processRegion(Region, VariantPrinter)` | L51-L60 | yes | Creates scope, runs `splicingPipeline`, joins |
| `createParallelMode()` | L66-L76 | yes | Creates anonymous parallel mode that submits `VardictWorker` per region |
| `VardictWorker(Region)` | L82-L97 | yes | `Callable<OutputStream>`: creates stream, runs `processRegion`, returns stream |
| `printHeader()` | L105-L110 | yes | Prints 4-column header: Sample, Chr, Intron, Intron count |

---

## Method Analyses

### AbstractMode.AbstractMode(segments, referenceResource)

**Source**: AbstractMode.java:L30-L33  
**Purpose**: Store input data for the mode  
**Called By**: All mode subclass constructors via `super()`  

**Algorithm**:
1. Store `segments` parameter (List of List of Region) into field `this.segments`
2. Store `referenceResource` parameter into field `this.referenceResource`

**Mutable State**: Fields `segments` and `referenceResource` are set once.

**Parity Notes**: None — trivial constructor.

---

### AbstractMode.pipeline(initialDataScope, executor)

**Source**: AbstractMode.java:L42-L57  
**Purpose**: Execute the FULL variant-calling pipeline for one region  
**Called By**: `SimpleMode.processBamInPipeline()`, `SomaticMode.processBothBamsInPipeline()`, `AmpliconMode.notParallel()`, `AmpliconMode.createParallelMode()`, `SomaticPostProcessModule.combineAnalysis()` (re-entrant)  

**Algorithm**:
1. Create a `CompletableFuture` via `supplyAsync()` that runs `new SAMFileParser().process(initialDataScope)` on the given `executor`
   - When executor is `DirectThreadExecutor`, this runs synchronously on the calling thread
   - Returns `Scope<VariationData>` (SAMFileParser output)
2. Chain `.thenApply(new CigarParser(false)::process)` — CigarParser with `splice=false`
   - Input: `Scope<VariationData>` from SAMFileParser
   - Output: `Scope<VariationData>` with populated variation/softclip maps
3. Chain `.thenApply(new VariationRealigner()::process)`
   - Input: `Scope<VariationData>` from CigarParser
   - Output: `Scope<VariationData>` with realigned variants
4. Chain `.thenApply(new StructuralVariantsProcessor()::process)`
   - Input: `Scope<VariationData>` from VariationRealigner
   - Output: `Scope<VariationData>` with structural variants detected
5. Chain `.thenApply(new ToVarsBuilder()::process)`
   - Input: `Scope<VariationData>` from StructuralVariantsProcessor
   - Output: `Scope<AlignedVarsData>` — final aligned variants map
6. Chain `.exceptionally()` handler:
   - Calls `stopVardictWithException(initialDataScope.region, ex)`
   - Then throws `new RuntimeException(ex)` to propagate

**Stage Typing**:
```
Scope<InitialData> → SAMFileParser → Scope<VariationData>
    → CigarParser(false) → Scope<VariationData>
    → VariationRealigner → Scope<VariationData>
    → StructuralVariantsProcessor → Scope<VariationData>
    → ToVarsBuilder → Scope<AlignedVarsData>
```

**Parity-Critical Observations**:
- Each stage creates a **new instance** (`new CigarParser(false)`, `new VariationRealigner()`, etc.) — no shared mutable state between pipeline stages except through `Scope`
- The `CigarParser` constructor receives `false` meaning this is NOT a partial pipeline call (affects SV-related logic within CigarParser)
- **Re-entrant invocations**: `StructuralVariantsProcessor` and `VariationRealigner` call `getMode().partialPipeline()` (not `pipeline()`); `SomaticPostProcessModule` calls `getMode().pipeline()` for combine analysis. All re-entrant calls use `DirectThreadExecutor`.
- Exception handler calls `System.exit(1)` — in Rust this should be process abort or error propagation

---

### AbstractMode.partialPipeline(currentScope, executor)

**Source**: AbstractMode.java:L79-L83  
**Purpose**: Run a PARTIAL pipeline (SAMFileParser + CigarParser only) for re-entrant structural variant / realignment region extension  
**Called By**: `StructuralVariantsProcessor` (5 call sites: L385, L555, L733, L1575, L1718), `VariationRealigner` (4 call sites: L1269, L1481, L1872, L2037)  

**Algorithm**:
1. Create `CompletableFuture` via `supplyAsync()` running `new SAMFileParser().process(currentScope)` on given executor
2. Chain `.thenApply(new CigarParser(true)::process)` — CigarParser with `splice=true`
   - Returns `Scope<VariationData>`
   - No further stages (no realigner, no SV processor, no ToVars)

**Stage Typing**:
```
Scope<InitialData> → SAMFileParser → Scope<VariationData>
    → CigarParser(true) → Scope<VariationData>
```

**Parity-Critical Observations**:
- `CigarParser(true)` — the `true` flag means this is a partial/re-entrant call. In CigarParser, this flag (called `isPartialPipeline` or `splice`) controls whether certain SV-related features are skipped or modified.
- The `currentScope` passed in often has **pre-populated variation maps** from the parent pipeline — CigarParser appends to them rather than starting fresh.
- Always invoked with `DirectThreadExecutor` — always synchronous.
- No `.exceptionally()` handler — exceptions propagate up to the calling pipeline stage's handler.

---

### AbstractMode.splicingPipeline(currentScope, executor)

**Source**: AbstractMode.java:L93-L97  
**Purpose**: Run a splicing-only pipeline (SAMFileParser + CigarParser) for splice-count output  
**Called By**: `SplicingMode.processRegion()`  

**Algorithm**:
1. Create `CompletableFuture` via `supplyAsync()` running `new SAMFileParser().process(currentScope)` on given executor
2. Chain `.thenApply(new CigarParser(false)::process)` — CigarParser with `splice=false`
   - Returns `Scope<VariationData>`

**Parity-Critical Observations**:
- Identical to `partialPipeline` except CigarParser receives `false` instead of `true`
- The difference: `splicingPipeline` uses `CigarParser(false)` = full cigar parsing (not partial), same as main pipeline; `partialPipeline` uses `CigarParser(true)` = partial/re-entrant mode
- No post-processing stages follow — splice counts are side-effected during CigarParser via the `splice` Set in Scope

---

### AbstractMode.stopVardictWithException(region, ex)

**Source**: AbstractMode.java:L64-L69  
**Purpose**: Fatal error handler — prints region, prints stack trace, exits process  
**Called By**: All pipeline `.exceptionally()` handlers, `tryToGetReference()`  

**Algorithm**:
1. Print to stderr: `"Critical exception occurs on region: {chr}:{start}-{end}, program will be stopped."`
2. Call `ex.printStackTrace()` to stderr
3. Call `System.exit(1)` — JVM process termination

**Parity Notes**: In Rust, this maps to `eprintln!` + `std::process::exit(1)` or potentially error propagation with `Result`. The exact error message format must match if any tests compare stderr output.

---

### AbstractMode.parallel()

**Source**: AbstractMode.java:L101-L103  
**Purpose**: Entry point for multi-threaded mode  
**Called By**: `VarDictLauncher.start()` when `threads > 1`  

**Algorithm**:
1. Call `createParallelMode()` (abstract — each subclass provides its own)
2. Call `.process()` on the returned `AbstractParallelMode` instance

---

### AbstractMode.tryToGetReference(region)

**Source**: AbstractMode.java:L143-L150  
**Purpose**: Safely get reference sequence for a region  
**Called By**: `SimpleMode.processBamInPipeline()`, `SomaticMode.notParallel()`, `SomaticMode.createParallelMode()`, `AmpliconMode.notParallel()`, `AmpliconMode.createParallelMode()`  

**Algorithm**:
1. Create empty `Reference reference = new Reference()`
2. Try: `reference = referenceResource.getReference(region)`
3. Catch any Exception: call `stopVardictWithException(region, ex)`
4. Return `reference`

**Edge Cases**:
- If `getReference()` throws, `stopVardictWithException` calls `System.exit(1)`, so line 4 is never reached in the exception case
- The empty `new Reference()` is a defensive default but never actually returned on success (overwritten by line 2)

---

### AbstractParallelMode.process()

**Source**: AbstractMode.java:L114-L134  
**Purpose**: Consumer loop for parallel mode — takes results from blocking queue, prints them  
**Called By**: `AbstractMode.parallel()`  

**Algorithm**:
1. Submit `produceTasks()` as a lambda to the executor (runs on a thread pool thread)
   - `produceTasks()` submits all region work items and ends by putting `LAST_SIGNAL_FUTURE` sentinel
2. Enter infinite consumer loop on the calling thread:
   3. `toPrint.take()` — blocks until a `Future<OutputStream>` is available
   4. If the taken future `== LAST_SIGNAL_FUTURE`, break out of loop
   5. Create a `VariantPrinter` via `VariantPrinter.createPrinter(instance().printerTypeOut)`
   6. Call `variantPrinter.print(wrk.get())` — blocks until the future completes, then prints the output
7. After loop: call `executor.shutdown()`

**Threading Model**:
- Producer thread: one thread from pool runs `produceTasks()`, which submits `Future`s to `toPrint` queue
- Consumer thread: the **calling thread** (main thread) takes futures from queue and prints in order
- Queue capacity: `CAPACITY = 10` — back-pressure if producer is more than 10 regions ahead
- Ordering: output order = region submission order (queue is FIFO), preserving BED file order

**Parity-Critical**:
- Output order is PRESERVED even in parallel mode because the consumer processes in FIFO order
- `wrk.get()` may throw `ExecutionException` wrapping pipeline exceptions — caught at line L132 and rethrown as `RuntimeException`
- `LAST_SIGNAL_FUTURE` identity check (`==`) is safe because it's a `static final` singleton

---

### SimpleMode.SimpleMode(segments, referenceResource)

**Source**: SimpleMode.java:L32-L35  
**Purpose**: Construct simple mode and print header  
**Called By**: `VarDictLauncher.start()` when no BAM2 and not amplicon/splicing  

**Algorithm**:
1. Call `super(segments, referenceResource)` — stores segments and referenceResource
2. Call `printHeader()` — outputs TSV header to stdout if `conf.printHeader` is true

---

### SimpleMode.notParallel()

**Source**: SimpleMode.java:L41-L47  
**Purpose**: Sequential processing of all regions in single-sample mode  
**Called By**: `VarDictLauncher.start()` when `threads == 1`  

**Algorithm**:
1. Create one `VariantPrinter` via `VariantPrinter.createPrinter(instance().printerTypeOut)`
2. For each `list` in `segments` (outer loop — typically one per chromosome):
   3. For each `region` in `list` (inner loop — intervals from BED file):
      4. Call `processBamInPipeline(region, variantPrinter)`

**Parity Notes**:
- Single `VariantPrinter` shared across all regions — output goes directly to stdout
- Region processing order: outer segment order × inner region order — matches BED file order exactly

---

### SimpleMode.createParallelMode()

**Source**: SimpleMode.java:L53-L63  
**Purpose**: Create parallel worker producer that submits one `VardictWorker` per region  
**Called By**: `AbstractMode.parallel()`  

**Algorithm** (anonymous `AbstractParallelMode.produceTasks()`):
1. For each `list` in `segments`:
   2. For each `region` in `list`:
      3. Create `VardictWorker(region)` — a `Callable<OutputStream>`
      4. Submit to `executor.submit(worker)` → get `Future<OutputStream>`
      5. Put the future into `toPrint` blocking queue (may block if queue full at 10)
6. Put `LAST_SIGNAL_FUTURE` into `toPrint` — signals consumer to stop

---

### SimpleMode.VardictWorker.call()

**Source**: SimpleMode.java:L72-L83  
**Purpose**: Process one region in an isolated output stream for parallel mode  
**Called By**: Thread pool executor  

**Algorithm**:
1. Create `ByteArrayOutputStream baos`
2. Wrap in `PrintStream out`
3. Create `VariantPrinter` and set its output to `out` (not stdout)
4. Call `processBamInPipeline(region, variantPrinter)` — full pipeline for this region
5. Close `out`
6. Return `baos` as the `OutputStream`

**Parity Notes**: Each parallel worker gets its own `ByteArrayOutputStream`, so output from different regions never interleaves. The consumer thread in `AbstractParallelMode.process()` serializes printing.

---

### SimpleMode.processBamInPipeline(region, out)

**Source**: SimpleMode.java:L90-L104  
**Purpose**: Core method — runs full pipeline for one region in simple (single-sample) mode  
**Called By**: `notParallel()`, `VardictWorker.call()`  

**Algorithm**:
1. `Reference reference = tryToGetReference(region)` — load reference sequence for region
2. Create `Scope<InitialData> initialScope`:
   - `bam`: `instance().conf.bam.getBam1()` — first (only) BAM file path
   - `region`: the input region
   - `reference`: from step 1
   - `referenceResource`: the shared reference resource
   - `maxReadLength`: `0` (initial — will be computed during SAMFileParser)
   - `splice`: `new HashSet<>()` — fresh empty set for splice tracking
   - `out`: the variant printer
   - `data`: `new InitialData()` — empty initial payload
3. Call `pipeline(initialScope, new DirectThreadExecutor())` — synchronous full pipeline
   - Returns `CompletableFuture<Scope<AlignedVarsData>>`
4. Chain `.thenAccept(new SimplePostProcessModule(out))` — post-process and print variants
5. Chain `.exceptionally()` — calls `stopVardictWithException` on failure
6. Call `.join()` — blocks until complete (since `DirectThreadExecutor` is used, this is already complete)

**Scope Construction Details**:
- `maxReadLength = 0`: this is the initial value; SAMFileParser will compute the actual max read length during BAM parsing
- `splice = new HashSet<>()`: a fresh Set per region that accumulates splice junction information during CigarParser; only `SomaticMode` shares the splice set between bam1 and bam2 pipelines
- `bam.getBam1()`: the single BAM file for simple mode

---

### SimpleMode.printHeader()

**Source**: SimpleMode.java:L106-L118  
**Purpose**: Print TSV header for simple mode output  
**Called By**: Constructor  

**Algorithm**:
1. If `instance().conf.printHeader` is false, return (do nothing)
2. Build header string by tab-joining 36 columns:
   `Sample, Gene, Chr, Start, End, Ref, Alt, Depth, AltDepth, RefFwdReads, RefRevReads, AltFwdReads, AltRevReads, Genotype, AF, Bias, PMean, PStd, QMean, QStd, MQ, Sig_Noise, HiAF, ExtraAF, shift3, MSI, MSI_NT, NM, HiCnt, HiCov, 5pFlankSeq, 3pFlankSeq, Seg, VarType, Duprate, SV_info`
3. If `instance().conf.crisprCuttingSite != 0`, append tab + `"CRISPR"` (37th column)
4. Print header to stdout via `System.out.println(header)`

**Parity Notes**:
- Column count is 36 (or 37 with CRISPR) — must match exactly
- Header uses `Utils.join("\t", ...)` — verify tab delimiter

---

### SomaticMode.SomaticMode(segments, referenceResource)

**Source**: SomaticMode.java:L33-L36  
**Purpose**: Construct somatic mode and print header  
**Called By**: `VarDictLauncher.start()` when hasBam2 is true  

**Algorithm**:
1. Call `super(segments, referenceResource)`
2. Call `printHeader()`

---

### SomaticMode.notParallel()

**Source**: SomaticMode.java:L42-L51  
**Purpose**: Sequential processing in somatic (tumor/normal) mode  
**Called By**: `VarDictLauncher.start()` when `threads == 1`  

**Algorithm**:
1. Create one `VariantPrinter`
2. For each `list` in `segments`:
   3. For each `region` in `list`:
      4. Create `splice = new ConcurrentHashSet<>()` — **new per region**
      5. `Reference ref = tryToGetReference(region)` — load reference once per region
      6. Call `processBothBamsInPipeline(variantPrinter, region, splice, ref)`

**Parity-Critical**:
- `ConcurrentHashSet` is used (not `HashSet`) even in single-threaded mode — this is for consistency with parallel mode but the concurrent behavior doesn't matter when single-threaded
- `splice` set is shared between tumor and normal BAM processing for the same region — this is **intentional**: splice junctions discovered in tumor BAM are available when processing normal BAM
- Reference is loaded ONCE per region and shared between both BAM pipelines

---

### SomaticMode.createParallelMode()

**Source**: SomaticMode.java:L57-L70  
**Purpose**: Create parallel producer for somatic mode  

**Algorithm** (anonymous `produceTasks()`):
1. For each `list` in `segments`:
   2. For each `region` in `list`:
      3. Create `splice = new ConcurrentHashSet<>()`
      4. `Reference ref1 = tryToGetReference(region)` — ref loaded **on producer thread**
      5. Submit `SomdictWorker(region, splice, ref1)` to executor → get `Future<OutputStream>`
      6. Put future into `toPrint` queue
7. Put `LAST_SIGNAL_FUTURE` sentinel

**Parity Notes**: Reference loading happens on the producer thread (running in the pool), not the consumer thread. The reference is passed into the worker by value.

---

### SomaticMode.SomdictWorker.call()

**Source**: SomaticMode.java:L79-L95  
**Purpose**: Process one region's tumor+normal pipeline in parallel mode  
**Called By**: Thread pool executor  

**Algorithm**:
1. Create `ByteArrayOutputStream baos`
2. Wrap in `PrintStream out`
3. Create `VariantPrinter`, set its output to `out`
4. Call `processBothBamsInPipeline(variantPrinter, region, splice, ref)` — processes both BAMs for this region
5. Close `out`
6. Return `baos`

---

### SomaticMode.processBothBamsInPipeline(variantPrinter, region, splice, ref)

**Source**: SomaticMode.java:L97-L125  
**Purpose**: **Core somatic method** — runs full pipeline on tumor BAM, then on normal BAM, then combines  
**Called By**: `notParallel()`, `SomdictWorker.call()`  

**Algorithm**:
1. Create `Scope<InitialData> initialScope1`:
   - `bam`: `instance().conf.bam.getBam1()` — **tumor BAM**
   - `region`: input region
   - `reference`: `ref` (shared reference loaded once)
   - `referenceResource`: shared reference resource
   - `maxReadLength`: `0` (initial)
   - `splice`: the shared splice set (will be populated during tumor pipeline)
   - `out`: `variantPrinter`
   - `data`: `new InitialData()`
2. Call `pipeline(initialScope1, new DirectThreadExecutor())` — full pipeline on tumor BAM
   - Returns `CompletableFuture<Scope<AlignedVarsData>>` as `bam1VariationsFuture`
3. **`bam1Variations = bam1VariationsFuture.join()`** — block and wait for tumor pipeline to complete
   - This extracts the `Scope<AlignedVarsData>` result
4. Create `Scope<InitialData> initialScope2`:
   - `bam`: `instance().conf.bam.getBam2()` — **normal BAM**
   - `region`: same region
   - `reference`: same `ref`
   - `referenceResource`: same
   - **`maxReadLength`: `bam1Variations.maxReadLength`** — carries over from tumor pipeline
   - `splice`: same shared splice set (now populated from tumor)
   - `out`: same printer
   - `data`: `new InitialData()`
5. Call `pipeline(initialScope2, new DirectThreadExecutor())` — full pipeline on normal BAM
   - Returns `CompletableFuture<Scope<AlignedVarsData>>` as `bam2VariationFuture`
6. Chain `bam2VariationFuture.thenAcceptBoth(bam1VariationsFuture, new SomaticPostProcessModule(referenceResource, variantPrinter))`
   - `thenAcceptBoth` receives BOTH results: normal (bam2) as first arg, tumor (bam1) as second arg
   - `SomaticPostProcessModule` compares tumor vs normal variants
7. Chain `.exceptionally()` — `stopVardictWithException`
8. Call `.join()` — block until somatic comparison is done

**Parity-Critical Observations**:
- **Processing order**: tumor BAM first, then normal BAM — this is REQUIRED because:
  - `maxReadLength` from tumor is used to initialize normal pipeline scope
  - `splice` set from tumor is available to normal pipeline
- **`thenAcceptBoth` argument order**: `bam2VariationFuture.thenAcceptBoth(bam1VariationsFuture, handler)` — the `BiConsumer` receives `(bam2_result, bam1_result)`. In `SomaticPostProcessModule.accept()`, the first arg is normal (bam2), second is tumor (bam1). This ordering is critical for correct somatic comparison.
- **Reference sharing**: Same `Reference` object used for both tumor and normal — this is safe because `Reference` is read-only after creation.
- **`SomaticPostProcessModule` re-entrant pipeline call**: The somatic post-processor itself calls `getMode().pipeline()` (at SomaticPostProcessModule.java:L404) for combine analysis regions. This is another full pipeline invocation using `DirectThreadExecutor`, creating a recursive pipeline call chain.

---

### SomaticMode.printHeader()

**Source**: SomaticMode.java:L127-L140  
**Purpose**: Print somatic mode TSV header  
**Called By**: Constructor  

**Algorithm**:
1. If `conf.printHeader` is false, return
2. Build header by tab-joining columns. The somatic header has **duplicated column groups** for tumor and normal:
   ```
   Sample, Gene, Chr, Start, End, Ref, Alt,
   [Tumor]: Depth, AltDepth, RefFwdReads, RefRevReads, AltFwdReads, AltRevReads, Genotype, AF, Bias, PMean, PStd, QMean, QStd, MQ, Sig_Noise, HiAF, ExtraAF, NM,
   [Normal]: Depth, AltDepth, RefFwdReads, RefRevReads, AltFwdReads, AltRevReads, Genotype, AF, Bias, PMean, PStd, QMean, QStd, MQ, Sig_Noise, HiAF, ExtraAF, NM,
   shift3, MSI, MSI_NT, 5pFlankSeq, 3pFlankSeq, Seg, VarLabel, VarType, Duprate1, SV_info1, Duprate2, SV_info2
   ```
3. Print to stdout

**Parity Notes**: 
- Somatic header has no `HiCnt`, `HiCov` columns (simple mode has them)
- Somatic header has `VarLabel` (not in simple mode)
- Somatic header has `Duprate1`, `SV_info1`, `Duprate2`, `SV_info2` (paired)
- The two repeated groups (tumor + normal) contain identical column names — Rust must reproduce this exactly

---

### AmpliconMode.AmpliconMode(segments, referenceResource)

**Source**: AmpliconMode.java:L36-L39  
**Purpose**: Construct amplicon mode  
**Called By**: `VarDictLauncher.start()` when `ampliconBasedCalling != null` and not splicing/simple/somatic  

**Algorithm**:
1. Call `super(segments, referenceResource)`
2. Call `printHeader()`

---

### AmpliconMode.notParallel()

**Source**: AmpliconMode.java:L45-L78  
**Purpose**: Sequential amplicon processing — processes all amplicons in a segment, then collects and post-processes  
**Called By**: `VarDictLauncher.start()` when `threads == 1`  

**Algorithm**:
1. Create one `VariantPrinter`
2. For each `regions` in `segments` (outer loop — each element is a list of overlapping amplicon regions):
   3. Create `pos = new HashMap<Integer, List<Tuple2<Integer, Region>>>()` — position→amplicon mapping
   4. `ampliconNumber = 0`
   5. `currentRegion = regions.get(0)` — initialize to first region
   6. Create `splice = new HashSet<>()` — shared across all amplicons in this segment
   7. Create `vars = new ArrayList<Map<Integer, Vars>>()` — collects variant maps from each amplicon
   8. For each `region` in `regions` (inner loop — one per amplicon):
      9. `currentRegion = region` — track the last region (will be used for post-processing)
      10. For each position `p` from `region.insertStart` to `region.insertEnd` (inclusive):
          11. `list = pos.computeIfAbsent(p, k -> new ArrayList<>())` — get or create position entry
          12. `list.add(tuple(ampliconNumber, region))` — add (amplicon index, region) to position map
      13. Create `Scope<InitialData> initialScope`:
          - `bam`: `instance().conf.bam.getBam1()`
          - `region`: current amplicon region
          - `reference`: `tryToGetReference(region)` — fresh reference per amplicon
          - `referenceResource`: shared
          - `maxReadLength`: `0`
          - `splice`: shared splice set
          - `out`: `variantPrinter`
          - `data`: `new InitialData()`
      14. Call `pipeline(initialScope, new DirectThreadExecutor())` → synchronous full pipeline
      15. `data = pipeline.join().data` — extract `AlignedVarsData`
      16. `vars.add(data.alignedVariants)` — collect the variant map
      17. `ampliconNumber++`
   18. Call `new AmpliconPostProcessModule().process(currentRegion, vars, pos, splice, variantPrinter)`
       - `currentRegion` is the **last** region in the amplicon list
       - `vars` contains aligned variants from ALL amplicons
       - `pos` maps genomic positions to which amplicons cover them

**Parity-Critical Observations**:
- **`currentRegion` is the LAST region**: After the loop, `currentRegion` holds the last element of `regions`. This is passed to `AmpliconPostProcessModule`, which uses it for gene name and chromosome info.
- **Position mapping**: `pos` maps each position in the *insert* range (not the full amplicon range) to the amplicon indices that cover it. `insertStart`/`insertEnd` are the "inner" amplicon boundaries (columns 7 and 8 from the BED file, 0-based).
- **Splice set is shared**: All amplicons in a segment share the same `splice` set — junctions discovered in one amplicon are visible to others.
- **Reference is loaded per amplicon**: Unlike somatic mode (one ref per region), amplicon mode loads reference per individual amplicon region.
- **Order matters**: Amplicons are processed in order, and `ampliconNumber` is assigned sequentially (0, 1, 2, ...). The post-processor relies on this ordering.

---

### AmpliconMode.createParallelMode()

**Source**: AmpliconMode.java:L84-L131  
**Purpose**: Parallel amplicon processing  
**Called By**: `AbstractMode.parallel()`  

**Algorithm** (anonymous `produceTasks()`):
1. For each `regions` in `segments`:
   2. Create `pos = new HashMap<>()` — same position mapping as sequential mode
   3. `j = 0` (amplicon counter)
   4. `currentRegion = regions.get(0)` — initialize
   5. Create `splice = new ConcurrentHashSet<>()` — concurrent version for parallel mode
   6. Create `workers = new ArrayList<CompletableFuture<Scope<AlignedVarsData>>>()` — collects pipeline futures

   7. For each `region` in `regions`:
      8. `currentRegion = region`
      9. Build position map (same as sequential: insertStart to insertEnd → tuple(j, region))
      10. Create `VariantPrinter` (per amplicon — for isolation in parallel)
      11. `Reference reference = tryToGetReference(region)` — load reference
      12. Create `Scope<InitialData> initialScope`:
          - Same as sequential mode but with fresh VariantPrinter
      13. Call `pipeline(initialScope, executor)` — **uses thread pool executor** (truly parallel)
          - Each amplicon pipeline can run on a different thread
      14. `workers.add(pipeline)`
      15. `j++`

   16. Collect results: for each `future` in `workers`:
       17. `vars.add(future.join().data.alignedVariants)` — blocks until each completes, collects in order

   18. Create `CompletableFuture<OutputStream> processAmpliconOutput = CompletableFuture.supplyAsync(...)`:
       - Creates `ByteArrayOutputStream` + `PrintStream`
       - Creates new `VariantPrinter` redirected to the byte stream
       - Calls `new AmpliconPostProcessModule().process(lastRegion, vars, pos, splice, variantPrinter)`
       - Returns the byte stream
       - With `.exceptionally()` handler
   19. `toPrint.add(processAmpliconOutput)` — note: uses `.add()` not `.put()` (non-blocking, but queue should have capacity since futures were already collected)
20. `toPrint.put(LAST_SIGNAL_FUTURE)` sentinel

**Parity-Critical**:
- Amplicon pipelines run in PARALLEL on the thread pool, but results are collected in ORDER via the `workers` list.
- Post-processing happens as an async task submitted to the same executor, but waits for all amplicon pipelines first.
- Uses `toPrint.add()` (non-blocking) for the amplicon output future and `toPrint.put()` (blocking) for the sentinel — the `.add()` could theoretically fail if queue is full, but shouldn't because the pipeline futures are already joined.
- The `splice` set is `ConcurrentHashSet` because multiple amplicon pipelines may write to it simultaneously.

---

### AmpliconMode.printHeader()

**Source**: AmpliconMode.java:L133-L142  
**Purpose**: Print amplicon mode TSV header  
**Called By**: Constructor  

**Algorithm**:
1. If `conf.printHeader` is false, return
2. Build header — same as simple mode base columns PLUS 4 amplicon-specific columns:
   ```
   Sample, Gene, Chr, Start, End, Ref, Alt, Depth, AltDepth, RefFwdReads, RefRevReads, AltFwdReads, AltRevReads, Genotype, AF, Bias, PMean, PStd, QMean, QStd, MQ, Sig_Noise, HiAF, ExtraAF, shift3, MSI, MSI_NT, NM, HiCnt, HiCov, 5pFlankSeq, 3pFlankSeq, Seg, VarType, GoodVarCount, TotalVarCount, Nocov, Ampflag
   ```
3. Print to stdout

**Parity Notes**:
- 38 columns total (34 + GoodVarCount, TotalVarCount, Nocov, Ampflag)
- No `Duprate` or `SV_info` columns
- No CRISPR column
- Column 35-38 are amplicon-specific: `GoodVarCount`, `TotalVarCount`, `Nocov`, `Ampflag`

---

### SplicingMode.SplicingMode(segments, referenceResource)

**Source**: SplicingMode.java:L29-L32  
**Purpose**: Construct splicing mode  
**Called By**: `VarDictLauncher.start()` when `outputSplicing == true`  

**Algorithm**:
1. Call `super(segments, referenceResource)`
2. Call `printHeader()`

---

### SplicingMode.notParallel()

**Source**: SplicingMode.java:L38-L44  
**Purpose**: Sequential splicing mode — iterate all regions  
**Called By**: `VarDictLauncher.start()` when `threads == 1`  

**Algorithm**:
1. Create one `VariantPrinter`
2. For each `list` in `segments`:
   3. For each `region` in `list`:
      4. Call `processRegion(region, variantPrinter)`

---

### SplicingMode.processRegion(region, out)

**Source**: SplicingMode.java:L51-L60  
**Purpose**: Run splicing-only pipeline for one region  
**Called By**: `notParallel()`, `VardictWorker.call()`  

**Algorithm**:
1. `Reference reference = tryToGetReference(region)`
2. Create `Scope<InitialData> initialScope`:
   - `bam`: `instance().conf.bam.getBam1()`
   - `region`: input region
   - `reference`: from step 1
   - `referenceResource`: shared
   - `maxReadLength`: `0`
   - `splice`: `new HashSet<>()` — fresh per region
   - `out`: variant printer
   - `data`: `new InitialData()`
3. Call `splicingPipeline(initialScope, new DirectThreadExecutor())`
   - Returns `CompletableFuture<Scope<VariationData>>`
4. Call `.join()` — block until complete

**Parity Notes**:
- No post-processing module — splicing output is a side-effect of CigarParser writing to the splice Set and/or printing splice counts directly
- Pipeline returns `Scope<VariationData>` but it's not used after `.join()`

---

### SplicingMode.createParallelMode()

**Source**: SplicingMode.java:L66-L76  
**Purpose**: Create parallel producer for splicing mode  

**Algorithm** (anonymous `produceTasks()`):
1. For each `list` in `segments`:
   2. For each `region` in `list`:
      3. Submit `VardictWorker(region)` to executor
      4. Put future into `toPrint` queue
5. Put `LAST_SIGNAL_FUTURE` sentinel

---

### SplicingMode.VardictWorker.call()

**Source**: SplicingMode.java:L82-L97  
**Purpose**: Process one region's splicing pipeline in parallel mode  

**Algorithm**:
1. Create `ByteArrayOutputStream baos`
2. Wrap in `PrintStream out`
3. Create `VariantPrinter`, set output to `out`
4. Call `processRegion(region, variantPrinter)` — splicing pipeline
5. Close `out`
6. Return `baos`

---

### SplicingMode.printHeader()

**Source**: SplicingMode.java:L105-L110  
**Purpose**: Print splicing mode header  
**Called By**: Constructor  

**Algorithm**:
1. If `conf.printHeader` is false, return
2. Tab-join 4 columns: `Sample, Chr, Intron, Intron count`
3. Print to stdout

**Parity Notes**: Only 4 columns — smallest header of all modes.

---

## Pipeline Composition per Mode

### Full Pipeline (used by Simple, Somatic, Amplicon)

```
SAMFileParser.process(Scope<InitialData>)          → Scope<VariationData>
CigarParser(false).process(Scope<VariationData>)    → Scope<VariationData>
VariationRealigner.process(Scope<VariationData>)    → Scope<VariationData>
StructuralVariantsProcessor.process(Scope<VariationData>) → Scope<VariationData>
ToVarsBuilder.process(Scope<VariationData>)         → Scope<AlignedVarsData>
```

Post-Processing:
- **Simple**: `→ SimplePostProcessModule.accept(Scope<AlignedVarsData>)`
- **Somatic**: `→ SomaticPostProcessModule.accept(Scope<AlignedVarsData> normal, Scope<AlignedVarsData> tumor)`
- **Amplicon**: `→ AmpliconPostProcessModule.process(lastRegion, List<alignedVariants>, posMap, splice, printer)`

### Partial Pipeline (re-entrant, used by StructuralVariantsProcessor and VariationRealigner)

```
SAMFileParser.process(Scope<InitialData>)           → Scope<VariationData>
CigarParser(true).process(Scope<VariationData>)     → Scope<VariationData>
```

- Always invoked via `getMode().partialPipeline(scope, DirectThreadExecutor)`
- `CigarParser(true)` = partial mode — may skip or modify SV-related logic
- Pre-populated variation maps in scope are augmented, not replaced

### Splicing Pipeline (used by SplicingMode only)

```
SAMFileParser.process(Scope<InitialData>)           → Scope<VariationData>
CigarParser(false).process(Scope<VariationData>)    → Scope<VariationData>
```

- Same as partial pipeline but `CigarParser(false)` — full parsing mode
- No post-processing; splice info is side-effected during parsing

### Somatic Combine Pipeline (re-entrant from SomaticPostProcessModule)

```
Full pipeline (same as above) invoked via:
    getMode().pipeline(Scope<InitialData>, DirectThreadExecutor).join()
```

- Called from `SomaticPostProcessModule` line 404
- Scope has `bam = bam1 + ":" + bam2` (concatenated BAM paths)
- Region is expanded: `variant.startPosition - maxReadLength` to `variant.endPosition + maxReadLength`
- Purpose: re-analyze an extended region to get combined coverage for a specific variant

---

## Cross-Module Dependencies

### Outbound Calls (Modes → other modules)

| Target Module | Method | Called From |
|--------------|--------|------------|
| `SAMFileParser` | `.process()` | `pipeline()`, `partialPipeline()`, `splicingPipeline()` |
| `CigarParser` | `.process()` | `pipeline()`, `partialPipeline()`, `splicingPipeline()` |
| `VariationRealigner` | `.process()` | `pipeline()` |
| `StructuralVariantsProcessor` | `.process()` | `pipeline()` |
| `ToVarsBuilder` | `.process()` | `pipeline()` |
| `SimplePostProcessModule` | `.accept()` | `SimpleMode.processBamInPipeline()` |
| `SomaticPostProcessModule` | `.accept()` | `SomaticMode.processBothBamsInPipeline()` |
| `AmpliconPostProcessModule` | `.process()` | `AmpliconMode.notParallel()`, `AmpliconMode.createParallelMode()` |
| `VariantPrinter` | `.createPrinter()`, `.setOut()`, `.print()` | All modes |
| `ReferenceResource` | `.getReference()` | `tryToGetReference()` |
| `GlobalReadOnlyScope` | `.instance()` | All modes (config access) |
| `DirectThreadExecutor` | (constructor) | All pipeline invocations in sequential mode |
| `Utils.join()` | tab-joining | `printHeader()` in all modes |

### Inbound Callers (other modules → Modes)

| Caller Module | Method Called | Purpose |
|--------------|-------------|---------|
| `VarDictLauncher` | `new *Mode()`, `setMode()`, `.notParallel()`, `.parallel()` | Launch mode |
| `StructuralVariantsProcessor` | `getMode().partialPipeline()` | Re-entrant SV region extension (5 call sites) |
| `VariationRealigner` | `getMode().partialPipeline()` | Re-entrant realignment region extension (4 call sites) |
| `SomaticPostProcessModule` | `getMode().pipeline()` | Re-entrant combine analysis (1 call site) |

---

## Known Parity Traps

1. **CigarParser `splice` flag semantics**: `pipeline()` creates `CigarParser(false)`, `partialPipeline()` creates `CigarParser(true)`, `splicingPipeline()` creates `CigarParser(false)`. This `boolean` flag is the ONLY difference between the three pipelines in terms of CigarParser behavior. In the Java code, the constructor parameter is named differently in different contexts ("splice", "isPartialPipeline"). Rust must match the exact flag value for each pipeline type.

2. **Somatic BAM processing order is tumor-first, then normal**: `processBothBamsInPipeline()` runs tumor (`bam1`) pipeline first, joins it, then runs normal (`bam2`) pipeline with tumor's `maxReadLength`. Reversing this order would produce different results because the splice set and maxReadLength would differ.

3. **`thenAcceptBoth` argument order in SomaticMode**: `bam2VariationFuture.thenAcceptBoth(bam1VariationsFuture, handler)` — the BiConsumer receives `(bam2_result, bam1_result)` = `(normal, tumor)`. `SomaticPostProcessModule` expects this exact order. If Rust inverts the arguments, tumor/normal columns will be swapped.

4. **Amplicon `currentRegion` is the LAST region in the list**: After iterating all amplicon regions, `currentRegion` holds the last entry. This is passed to `AmpliconPostProcessModule.process()` as the "representative" region. The last region's gene name and chromosome are used for output.

5. **Amplicon position mapping uses `insertStart`/`insertEnd` (not `start`/`end`)**: The `pos` map ranges over `region.insertStart..=region.insertEnd`, which are the inner amplicon boundaries (columns 7-8 from BED). Using `region.start`/`region.end` instead would produce wrong position→amplicon assignments.

6. **Splice set scope differs by mode**:
   - **SimpleMode**: `new HashSet<>()` — fresh per region, not shared
   - **SomaticMode**: `new ConcurrentHashSet<>()` — shared between tumor and normal for same region
   - **AmpliconMode (sequential)**: `new HashSet<>()` — shared across ALL amplicons in a segment
   - **AmpliconMode (parallel)**: `new ConcurrentHashSet<>()` — shared across all amplicons in a segment
   - **SplicingMode**: `new HashSet<>()` — fresh per region
   
   Using the wrong splice set scope will cause parity mismatches in modes where splice information crosses BAM boundaries or amplicon boundaries.

7. **Re-entrant `partialPipeline` always uses `DirectThreadExecutor`**: Even when the outer mode is parallel, all re-entrant calls from `StructuralVariantsProcessor` and `VariationRealigner` use `DirectThreadExecutor` (synchronous). Rust must not accidentally dispatch these sub-pipelines to a thread pool.

8. **Re-entrant `pipeline` from `SomaticPostProcessModule` uses concatenated BAM path**: The combine analysis scope has `bam = bam1 + ":" + bam2` (colon-separated). This combined BAM path is passed to `SAMFileParser`, which must handle it (reading from both tumor and normal).

9. **`maxReadLength` propagation in somatic mode**: Tumor pipeline starts with `maxReadLength=0` and computes the actual value during SAMFileParser. This computed value is passed to the normal pipeline scope. In simple/amplicon modes, `maxReadLength` starts at 0 and is not propagated to other pipeline runs.

10. **Parallel mode output ordering is preserved**: `AbstractParallelMode.process()` consumes `Future`s from a FIFO `LinkedBlockingQueue`, which preserves submission order. Even though regions run in parallel, output rows appear in BED-file order. Rust must preserve this ordering guarantee.

11. **Producer submits to queue capacity 10 with back-pressure**: `LinkedBlockingQueue(CAPACITY=10)` means the producer blocks when 10 futures are already queued. This limits memory usage (at most 10 buffered region outputs). Rust should implement similar back-pressure to avoid unbounded memory for large BED files.

12. **Exception handling terminates the JVM**: `stopVardictWithException` calls `System.exit(1)`. In parallel mode, this kills all threads. Rust must decide whether to mimic this (process abort) or use error propagation.

13. **`printHeader()` called in constructor**: Every mode's constructor calls `printHeader()` before any pipeline work begins. The header is printed to stdout via `System.out.println`. In Rust, this must happen during mode construction, not after pipeline completion.

14. **Simple mode header conditionally adds CRISPR column**: When `conf.crisprCuttingSite != 0`, a 37th column `"CRISPR"` is appended. No other mode has this conditional column. Rust must check this condition and append the column.

15. **AmpliconMode parallel uses `.add()` vs `.put()`**: In `createParallelMode()`, the amplicon output future is added with `toPrint.add()` (non-blocking, throws if full) while the sentinel uses `toPrint.put()` (blocking). Since all amplicon pipeline futures are already `.join()`ed before adding, the queue should have capacity. But this is a subtle API difference from other modes that use `.put()` throughout.

16. **Somatic header has no `HiCnt`/`HiCov` columns**: Simple mode header includes `HiCnt` and `HiCov` between `NM` and `5pFlankSeq`. Somatic mode omits them. Amplicon mode includes them. Column counts differ:
    - Simple: 36 (+ optional CRISPR = 37)
    - Somatic: 52 (fixed)
    - Amplicon: 38 (fixed)
    - Splicing: 4 (fixed)

17. **`Reference` object lifecycle differs by mode**: 
    - SimpleMode: loaded per region, used by one pipeline  
    - SomaticMode: loaded per region, shared by tumor+normal pipelines  
    - AmpliconMode: loaded per amplicon region (not per segment)
    - This affects memory usage but not correctness unless reference loading has side-effects

18. **`GlobalReadOnlyScope.setMode()` is set ONCE**: The mode singleton is set by `VarDictLauncher.start()` after construction. It's accessed via `getMode()` by `StructuralVariantsProcessor`, `VariationRealigner`, and `SomaticPostProcessModule`. In Rust, this could be a shared reference or a global — but it must be available during re-entrant pipeline calls.

19. **AmpliconMode parallel collects amplicon results BEFORE post-processing**: In parallel mode, ALL amplicon pipelines are launched, then their results are collected (joined) in order, then post-processing runs. This differs from sequential mode where each amplicon is processed and collected one-by-one. The functional result should be identical since amplicon pipelines are independent (sharing only the splice set).

20. **SplicingMode pipeline returns `Scope<VariationData>`, not `Scope<AlignedVarsData>`**: The splicing pipeline skips `VariationRealigner`, `StructuralVariantsProcessor`, and `ToVarsBuilder`. The returned scope type is `VariationData` (not `AlignedVarsData`), and the result is discarded after `.join()`. All splicing output is side-effected during CigarParser.
