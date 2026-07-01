# Known Issues

## `inic()` dispatch mismatches

These issues are recorded only; they are not corrected in this change.

| Requested algorithm | File and line | Observed dispatch | Expected dispatch |
| --- | --- | --- | --- |
| `1Derberaprox_v2` | `src/MaGICS/init.c:253-255` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Derberaprox_v2+` | `src/MaGICS/init.c:260-263` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Derberaprox_v2-` | `src/MaGICS/init.c:268-271` | `propagate_1D_erberaprox` | `propagate_1D_erberaprox_v2` |
| `1Dklepikov_v2` | `src/MaGICS/init.c:299-302` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |
| `1Dklepikov_v2+` | `src/MaGICS/init.c:306-309` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |
| `1Dklepikov_v2-` | `src/MaGICS/init.c:314-317` | `propagate_1D_erberaprox` | `propagate_1D_klepikov_v2` |

The default branch for no arguments does select `propagate_1D_klepikov_v2` at `src/MaGICS/init.c:221-225`, so the mismatch affects explicit `*_v2` arguments.
