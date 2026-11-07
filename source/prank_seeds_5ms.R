# Run PRANK with many random seeds on each directly aligned set of the 5 ms song sequences.
# PRANK gives different alignments with different seeds, so one alignment per set is chosen afterwards.
#
# Run from the project root:
#   Rscript source/prank_seeds_5ms.R <n seeds> <cores>

args <- commandArgs(TRUE)
n_seeds <- if (length(args) > 0) as.integer(args[1]) else 100
cores <- if (length(args) > 1) as.integer(args[2]) else 24

seq_dir <- "./data/processed/song_sequences_5ms/song_type_sequence"
out_dir <- "./data/processed/mafft_alignments_5ms/prank_seeds"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
work <- tempfile("prank_seeds_")
dir.create(work)

read_fa <- function(file) {
  a <- readLines(file)
  starts <- grep("^>", a)
  s <- mapply(function(i, j) paste(a[(i + 1):j], collapse = ""), starts, c(starts[-1] - 1, length(a)))
  setNames(as.vector(s), sub("^>", "", a[starts]))
}

fa <- list.files(seq_dir, pattern = "^song_seqs_(early_new|early_old|no\\.fossils)_[A-Z0-9]+_prank\\.fa$")
jobs <- expand.grid(file = fa, seed = seq_len(n_seeds), stringsAsFactors = FALSE)

runs_file <- file.path(out_dir, "prank_seed_runs.RDS")
t0 <- Sys.time()
if (!file.exists(runs_file)) {
o <- parallel::mclapply(seq_len(nrow(jobs)), mc.cores = cores, mc.preschedule = FALSE, function(i) {
  base <- file.path(work, paste0("job", i))
  log <- paste0(base, ".log")
  system(paste0("prank -d=", file.path(seq_dir, jobs$file[i]), " -o=", base,
                " -iterate=100 -kappa=2 -seed=", jobs$seed[i], " > ", log, " 2>&1"))
  aln <- read_fa(paste0(base, ".best.fas"))
  lines <- unlist(strsplit(readLines(log, warn = FALSE), "\r"))
  scores <- as.numeric(sub("^.*Alignment score: *", "", grep("Alignment score:", lines, value = TRUE)))
  unlink(paste0(base, c(".best.fas", ".log")))
  list(file = jobs$file[i], seed = jobs$seed[i], scores = scores, aln = aln[order(names(aln))])
})
bad <- vapply(o, inherits, TRUE, "try-error")
if (any(bad)) stop(sum(bad), " runs failed\n", paste(unique(o[bad]), collapse = "\n"))

saveRDS(o, runs_file)
}
unlink(work, recursive = TRUE)
o <- readRDS(runs_file)
cat("runs:", length(o), "\n")

# Selection: the most central alignment, the one that shares the most aligned residue pairs with the
# alignments from all other seeds (lowest seed if tied). PRANK's own score is kept for comparison: it is
# a penalty (each run outputs its lowest-scoring iteration)

# aligned residue pairs as integer keys: residue r of song i aligned with residue s of song j
pairs_of <- function(a) {
  m <- do.call(rbind, strsplit(a, ""))
  r <- t(apply(m != "-", 1, cumsum))
  r[m == "-"] <- 0
  unlist(lapply(seq_len(ncol(m)), function(cl) {
    w <- which(r[, cl] > 0)
    if (length(w) < 2) return(NULL)
    p <- combn(w, 2)
    (p[1, ] - 1) * 1e6 + (p[2, ] - 1) * 1e4 + (r[p[1, ], cl] - 1) * 100 + (r[p[2, ], cl] - 1)
  }))
}

sets <- split(o, sapply(o, `[[`, "file"))
selection <- parallel::mclapply(sets, mc.cores = min(cores, length(sets)), mc.preschedule = FALSE, function(x) {
  seeds <- sapply(x, `[[`, "seed")
  x <- x[order(seeds)]
  seeds <- sort(seeds)
  score <- sapply(x, function(z) min(z$scores))
  gaps <- sapply(x, function(z) nchar(z$aln[1]) - max(nchar(gsub("-", "", z$aln))))
  prop <- sapply(x, function(z) sum(nchar(gsub("[^-]", "", z$aln))) / sum(nchar(z$aln)))
  p <- lapply(x, function(z) pairs_of(z$aln))
  n <- length(p)
  sim <- diag(n)
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    k <- length(intersect(p[[i]], p[[j]]))
    sim[i, j] <- sim[j, i] <- k / (length(p[[i]]) + length(p[[j]]) - k)
  }
  central <- function(idx) idx[which.max(rowMeans(sim[idx, idx, drop = FALSE]))]
  lowest <- function(idx) idx[order(score[idx], gaps[idx], idx)[1]]
  half1 <- seq_len(n %/% 2)
  half2 <- setdiff(seq_len(n), half1)
  b <- central(seq_len(n))
  label <- sub("^song_seqs_(.*)_prank\\.fa$", "\\1", x[[1]]$file)
  aln <- x[[b]]$aln
  writeLines(as.vector(rbind(paste0(">", names(aln)), aln)), file.path(out_dir, paste0("selected_", label, ".fa")))
  data.frame(set = label, seeds = n, distinct_alignments = length(unique(sapply(x, function(z) paste(z$aln, collapse = "|")))),
             mean_similarity = mean(sim[upper.tri(sim)]), selected_seed = seeds[b], gaps = gaps[b], gaps_min = min(gaps),
             gaps_median = median(gaps), gaps_max = max(gaps), gap_prop = prop[b], score = score[b], lowest_score = min(score),
             centrality = mean(sim[b, ]), split_half_central = sim[central(half1), central(half2)],
             split_half_lowest_score = sim[lowest(half1), lowest(half2)], stringsAsFactors = FALSE)
})
selection <- do.call(rbind, selection)
rownames(selection) <- NULL
write.csv(selection, file.path(out_dir, "prank_selection.csv"), row.names = FALSE)
print(selection[, c("set", "distinct_alignments", "selected_seed", "gaps", "gaps_median", "mean_similarity", "split_half_central", "split_half_lowest_score")], digits = 2)
print(Sys.time() - t0)
