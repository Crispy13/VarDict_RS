# StructuralVariantsProcessor

**Source**: `src/mods/structural_variants_processor.rs`
**LOC**: 4,702
**Java counterpart**: `modules/StructuralVariantsProcessor.java`
**Status**: partial

## Overview
Owns structural-variant discovery after realignment. The hot `findsv()` path scans unused 5' and 3' soft clips, attempts forward DEL/DUP-style matches first, then falls back to reverse-complement INV discovery when forward matching fails.

## Public API
| Function/Method | Purpose |
|----------------|---------|
| `find_svs(...)` | Entry point for SV detection across the region |
| `find_svs_del_candidates(...)` | Java `findsv()` port for 5'/3' soft-clip DEL/DUP/INV candidate creation |
| `find_del_disc(...)` | Processes discordant-pair deletion clusters |

## Java Correspondence
| Rust | Java | Notes |
|------|------|-------|
| `find_svs_del_candidates()` | `StructuralVariantsProcessor.findsv()` | Rust keeps the direct forward match historical-window blind to avoid chr13/054 false positives |
| `check_pairs()` | `findsv()` DEL-pairs overlap scan | Marks matched clusters used, so synthetic probes must not call it directly |
| `is_softp2sv_first_used()` | `SOFTP2SV{p}->[0].used` checks | Rust reconstructs the Java bucket top-entry lazily from cluster vectors |

## Known Parity Traps
- `findsv()` forward matches cannot globally enable historical windows in Rust. Doing so reintroduces chr13/054 extra INV candidates.
- Java's mutable reference seed map grows as remote windows are loaded. Rust emulates this only on the synthetic pre-INV forward probe by temporarily loading the reverse partner span into history, then restoring the exact pre-probe active window and truncating helper-added history.
- The synthetic pre-INV pairs gate must be bucketed by `SoftClip.softp` to match Java `SOFTP2SV{softp}` behavior. A plain overlap scan across all DUP/DEL clusters can incorrectly retain spurious INVs by counting unrelated clusters.

## Divergences from Java
- Rust keeps a single active contiguous reference window plus historical window snapshots instead of Java's monotonically growing mutable `REF` hash. Parity-sensitive probes may need to materialize remote spans explicitly before seed lookups.
- Rust reconstructs several Java map views (`SOFTP2SV`, top-entry used state) on demand from the cluster vectors instead of persisting every intermediate bucket.

## Cross-Module Dependencies
- Calls: shared reference loading helpers, `Variant`/`SoftClip`, realigned variation maps, output SV counters.
- Called by: `variant_realigner.rs`, `vardict_pipeline.rs`.