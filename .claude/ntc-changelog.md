# TSScall — revision changelog (for the NTC ecosystem manager)

<!-- THIS FILE IS OWNED BY THIS REPO. The ecosystem manager (ntc-ecosystem-context) READS it to
track pipeline revisions and keep its overseer context current; it does NOT overwrite it.
Append a new entry at the TOP whenever you make a notable change (a revision / version bump /
results-affecting fix), so the manager has visibility. Keep entries to one or two lines. -->

Entry format (newest first):
```
## <YYYY-MM-DD> · <version/tag or claude/branch> · <NTC-### if any>
- <one-line summary>; affects results: <yes|no>
```

<!-- entries below, newest first -->

## 2026-11-26 · claude/vibrant-clarke-6uqmti
- PERFORMANCE.md: added real-data 3-way benchmark (Armstrong005 hg38) — orig==refac byte-identical (8/8 completed cells), 10–184× speedup, raw-original time-bound (~46 GB, >24 h DNF); refac/improved do raw in minutes. affects results: no (doc)
