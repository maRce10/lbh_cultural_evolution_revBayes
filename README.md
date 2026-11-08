Understanding birdsong cultural evolution through the fossilized
birth-death process
================

[Marcelo Araya-Salas](https://marce10.github.io/), [Beatriz
Willink](https://scholar.google.com/citations?user=0a8k9T8AAAAJ&hl=es&oi=ao)
and Alejandro Rico-Guevara

<!-- Description  -->

Creating long-billed hermit song cultural evolution trees and analyzing
evolutionary patterns

**Website:** [Modelling birdsong cultural
evolution](https://marce10.github.io/lbh_cultural_evolution_revBayes/)

*Updated on 2026-11-07 11:56:47.240031*

<!-- README.md is generated from README.Rmd. Please edit that file -->

## Table of contents

- [General info](#general-info)
- [Song segmenting](#song-segmenting)
- [Song alignment](#song-alignment)
- [Next steps](#next-steps)
- [Analysis](#analysis)
- [In progress](#in-progress)
- [To-do list](#to-do-list)
- [Website](#website)
- [Contact](#contact)

## General info

Code used for creating long-billed hermit song cultural evolution trees
and analyzing evolutionary patterns.

## Song segmenting

Songs are now coded as sequences of 5 ms segments with 5 possible
states, in four steps. The [full
report](https://marce10.github.io/lbh_cultural_evolution_revBayes/source/song_segment_coding_5ms.html)
is built from `source/song_segment_coding_5ms.qmd`. Measurements use
dominant frequency contours extracted with warbleR’s `freq_ts()` and
manually tailored where needed.

1.  **Manual broad sections**: each song is divided into pure tone and
    trill sections. Consecutive pure tone sections are merged, and
    trills that differ discretely in rate stay as separate sections
2.  **Classify trills**: each whole trill section is fast (F, period \<
    13 ms) or slow (S, period ≥ 13 ms)
3.  **Split into 5 ms segments**: the whole song is split, and trill
    segments take the state of their section
4.  **Classify pure tones**: each tone segment is descending (D, slope
    \< -14 kHz/s), flat (P) or ascending (U, slope ≥ 14 kHz/s), with the
    slope measured in a 15 ms window

Output files are in `data/processed`:

| Output | Five leks used in the paper (197 recordings, 104 song types) | All leks (322 recordings, 187 song types) |
|:---|:---|:---|
| Segments | `song_5ms_segments_by_recording.csv` | `song_5ms_segments_by_recording_all_leks.csv` |
| Sequence per song type and year | `song_5ms_sequences_by_taxon.csv` | `song_5ms_sequences_by_taxon_all_leks.csv` |
| One sequence per song type | `song_5ms_sequences_by_song_type.csv` | `song_5ms_sequences_by_song_type_all_leks.csv` |
| FASTA of all recordings | `song_sequences_5ms/all_recordings.fa` | `song_sequences_5ms/all_recordings_all_leks.fa` |
| FASTA files for alignment | `song_sequences_5ms/song_type_sequence/` and `song_sequences_5ms/year_specific_sequence/` | not available |

## Song alignment

The songs of the five leks in the paper are aligned using one sequence
per song type (`song_sequences_5ms/song_type_sequence/`), as in the
first version of the study. The [full
report](https://marce10.github.io/lbh_cultural_evolution_revBayes/source/song_alignment_5ms.html)
is built from `source/song_alignment_5ms.qmd`.

There are three alignment strategies:

| Strategy | File name | Method |
|:---|:---|:---|
| MAFFT-optimal | `optimal` | MAFFT L-INS-i with optimized scores: -2 (different category), -0.75 (same category), 2 (same state) |
| MAFFT-agnostic | `all.equal` | MAFFT L-INS-i with the same score for any mismatch: 1 (mismatch), 1.5 (same state) |
| PRANK | `prank` | The most central alignment among 100 seeds |

Each song type is aligned once, and the datasets are built from that
alignment, so `all` and `bounds` have identical aligned rows:

| Dataset | Songs included | Leks | Nexus files (× 3 strategies) | RevBayes models (× 2 clocks) |
|:---|:---|:---|---:|---:|
| `all_new` | All years, recordings after 2000 | BR1, CCE, HC1, SUR, TR1 | 15 | 30 |
| `bounds_new` | First and last year of each song type, recordings after 2000 | BR1, CCE, HC1, SUR, TR1 | 15 | 30 |
| `all_old` | All years, including earlier recordings | CCE, HC1, SUR | 9 | 18 |
| `bounds_old` | First and last year of each song type, including earlier recordings | CCE, HC1, SUR | 9 | 18 |
| `no.fossils` | Song types in the last year of sampling | BR1, SUR, TR1 | 9 | 18 |
| **Total** |  |  | **57** | **114** |

`bounds` replaces `early` (only the first year of each song type) from
the first version of the study. `early` misses song types that persist
until the last year of sampling, which conflicts with assuming that all
songs in that year were sampled (rho = 1), so it is now used only as the
base alignment. Files named `early` are not analyzed.

Output files and scripts:

| Output | Location | Script |
|:---|:---|:---|
| Nexus files for RevBayes | `data/processed/nexus_5ms/` | `source/write_nexus_5ms.R` |
| Alignments, scoring matrices and results | `data/processed/mafft_alignments_5ms/` | `source/align_song_types_5ms.R` |
| PRANK runs and selection | `data/processed/mafft_alignments_5ms/prank_seeds/` | `source/prank_seeds_5ms.R` |
| MAFFT score optimization | `data/processed/mafft_optimization_5ms/` | `source/optimize_mafft_costs_5ms.R` |

## Next steps

The analyses below were run on the earlier sequences
(`segments_by_song_type.csv`, 20 segments, 6 states) and **need to be
redone** with the new 5 ms, 5-state sequences, in this order.

**Sequences and alignments**

- [x] Choose the FASTA version for the alignments: `song_type_sequence`,
  as in the first version of the study
- [x] Redo the MAFFT cost optimization for the 5 states
- [x] Rerun the MAFFT and PRANK alignments, and check gaps and
  convergence
- [x] Replace `early` with `bounds` and write the nexus files
  (`data/processed/nexus_5ms/`)

**RevBayes**

- [ ] Point the `.Rev` scripts to `data/processed/nexus_5ms/`, with the
  new file names (e.g. `BR1_optimal_bounds_new_alignment.nex`)
- [ ] Change the GTR model from 6 to 5 states, which means 5 stationary
  frequencies and 10 exchangeability rates (`source/model_GTRG.Rev`,
  `source/no_fossils/model_GTRG_BD.Rev`)
- [ ] Confirm that rho = 1 and the treatment of gaps still fit the new
  data (about half of each alignment is gaps)
- [ ] Rerun all 114 models (`source/*.Rev`,
  `source/trees_with_revbayes.Rmd`), including those with no fossils
- [ ] Rerun the posterior predictions (P3; `source/Write_P3.Rmd`,
  `source/run_P3_*.sh`)

**Results**

- [ ] Redo the MCMC, ML and P3 diagnostics
  (`source/Model diagnostics.Rmd`, `source/rwty_diagnostics.Rmd`)
- [ ] Redo the tree topology diagnostics
  (`source/tree_topology_diagnostics_oct_2024.qmd`)
- [ ] Redo the descriptive statistics (`source/descriptive_stats.Rmd`)
- [ ] Redo the regression models (`source/regression_models.qmd`,
  `source/song_sequence_similarity.qmd`)
- [ ] Assess the effect of the redundancy of 5 ms sequences
  (neighbouring segments mostly repeat the same state) on the trees
- [ ] Check which results depend on the alignment strategy, and
  interpret results found only with PRANK with care, especially in BR1,
  SUR and the `old` datasets, where PRANK alignments vary among seeds

**Manuscript**

- [ ] Rewrite the methods on sound classification and song coding (now 5
  ms segments and 5 states instead of 20 segments and 6 categories),
  including Fig. 2
- [ ] Update the methods on alignment: new MAFFT scores and how they
  were selected, PRANK selection among seeds (both rules were set after
  inspecting the results), later years taking the aligned row of their
  song type, MAFFT convergence, and `bounds` instead of `early`
- [ ] Correct the leks in the analyses with no fossils: the text says “3
  or more distinct songs” in the last year, but the code uses more than
  3, which excludes HC1 (exactly 3)
- [ ] Replace the repeatability results with the agreement between
  copies of the same song type, and update the supplementary text
  (`source/segment_classification_repeatability.qmd`)
- [ ] Redo the figures (`manuscript/Final_Figures.Rmd`) and update all
  results
- [ ] Revise the discussion, including the assumption on discretizing
  songs of variable duration

**Website**

- [x] Add the song alignment section
- [ ] Add the RevBayes section once it is redone

## Analysis

Items marked *to be redone* use the earlier sequences (see [Next
steps](#next-steps)).

- [Song
  segmenting](https://marce10.github.io/lbh_cultural_evolution_revBayes/source/song_segment_coding_5ms.html)
- [Song
  alignment](https://marce10.github.io/lbh_cultural_evolution_revBayes/source/song_alignment_5ms.html)
- Revbayes models done (*to be redone*)
- Making trees with RevBayes (*to be redone*)
- Run revbayes posterior predictions (P3) (*to be redone*)
- Summarize results revBayes trees (*to be redone*)
- [MCMC, ML and P3
  diagnostics](https://rpubs.com/marcelo-araya-salas/722260) (*to be
  redone*)
- [diagnostic
  plots](lbh_cultural_evolution/tree/master/output/MCMC_diagnostic_plots)
  (*to be redone*)
- [tree topology
  diagnostics](https://rpubs.com/marcelo-araya-salas/709692) (*to be
  redone*)

## In progress

- Tree topology diagnostics

## To-do list

- New github without output (MAS)
- Add recording methods (MAS)
- Add tree congruence results to github (MAS)
- Add tree congruence to methods (MAS)
- Add Alejo comments to ms
- Run RevBayes with no fossils
- Check aligment convergence

## To-do list 3-12-2024

Marcelo:

- combinar graficos de distribuciones posteriores (decidir si dejar
  reloj)
- hacer figura 1 (poner espectrograma con panel con secuencia de sonidos
  sobre el espectrograma)
- materiales suplementarios con ejemplos de codigo para alineamientos y
  modelos
- hacer un nuevo repositorio (leer el manuscrito y asegurarse q los
  codigos de cada cosa estan representados)
- Agregar texto suplementario sobre como se midio repetibilidad y
  agregar valor de repitilibidad al texto principal
- meter en materiales suppl diagnosticos de modelos de reconstruccion
  (<https://rpubs.com/marcelo-araya-salas/722260>)
- Rhat y prsf son lo mismo, usar solo uno

## Website

The website is built with Quarto (`_quarto.yml`, `index.qmd`) and
deployed to GitHub Pages by `.github/workflows/publish.yml` on each push
to `master`. Code results are read from `_freeze/`, so after editing a
report run `quarto render` and commit `_freeze/`.

## Status

Project is: *in progress*

## Contact

Created by [Marcelo Araya-Salas](https://marceloarayasalas.weebly.com/)
