# Archive

Scripts, data and results of the first version of the study (songs coded as 20 segments with 6 states), and of
exploratory analyses that are not part of the current pipeline. The folder structure mirrors the project root
(`archive/source`, `archive/data`, `archive/output`; the current `scripts/` folder was called `source/`), so a file that was at `source/x.Rmd` is now at
`archive/source/x.Rmd`. Paths inside these files were not updated: to run one, move it (and the files it reads)
back to its original location.

Naming in the first version, and the equivalent in the current analyses:

| First version | Current |
|:---|:---|
| `old` (historical records including those before 2000) | `full` |
| `new` (records after 2000) | `recent` |
| `no.fossils` (only the song types of the last year of sampling) | `none` |
| `early` (first year of each song type) | `bounds` (first and last year); `early` is only the base alignment |
| `data/processed/nexus/`, `data/processed/fossils/` (now in `archive/`) | `data/processed/nexus/`, `data/processed/taxa/` |
| `all.equal` (MAFFT alignment with the same score for any mismatch) | `agnostic` |
| `source/` | `scripts/` |
