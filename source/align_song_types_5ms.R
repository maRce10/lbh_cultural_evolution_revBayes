# Align the 5 ms song sequences of the five paper leks.
#
# Follows source/aligning_song_types.Rmd, using one sequence per song type
# (song_sequences_5ms/song_type_sequence) and the same 57 alignments:
#   strategies: optimal (MAFFT, optimized scores), all.equal (MAFFT, same score for any mismatch), prank
#   early_new, early_old, no.fossils: aligned directly
#   all_new, all_old: the remaining years are added to the early alignment. All years of a song type have
#   the same sequence, so each year takes the aligned row of its song type (mafft --add realigns each copy
#   on its own and can place identical songs in different columns or drop part of them)
#   bounds_new, bounds_old: only the first and last year of each song type, taken from the all alignment.
#   The analyses compare all vs. bounds, which share the same aligned rows; early is only the base alignment
#
# Run from the project root:
#   Rscript source/align_song_types_5ms.R <cores>

args <- commandArgs(TRUE)
cores <- if (length(args) > 0) as.integer(args[1]) else 16

# same MAFFT version used to optimize the scores
mafft <- path.expand("~/mambaforge/envs/mafft/bin/mafft")
if (!file.exists(mafft)) stop("mafft not found")

seq_dir <- "./data/processed/song_sequences_5ms/song_type_sequence"
out_dir <- "./data/processed/mafft_alignments_5ms"
work <- tempfile("align_5ms_")
dir.create(work)
dir.create(out_dir, showWarnings = FALSE)

# scores for the all.equal strategy, as in the previous study (the optimal matrix comes from
# source/song_alignment_5ms.qmd)
hex <- c("0x64", "0x66", "0x70", "0x73", "0x75")
pairs <- expand.grid(a = hex, b = hex, stringsAsFactors = FALSE)
writeLines(sprintf("%s %s %.2f", pairs$a, pairs$b, ifelse(pairs$a == pairs$b, 1.5, 1)),
           file.path(out_dir, "all.equal_cost_matrix.txt"))

read_fa <- function(file) {
  a <- readLines(file)
  starts <- grep("^>", a)
  s <- mapply(function(i, j) if (j > i) paste(a[(i + 1):j], collapse = "") else "", starts, c(starts[-1] - 1, length(a)))
  setNames(as.vector(s), sub("^>", "", a[starts]))
}

# one sequence per line, whatever the aligner does
write_fa <- function(s, file) writeLines(as.vector(rbind(paste0(">", names(s)), s)), file)

align_direct <- function(job) {
  output <- file.path(out_dir, paste0(job$lek, "_", job$strategy, "_", job$set, "_alignment.fa"))
  if (job$strategy == "prank") {
    # selected among 100 seeds in source/prank_seeds_5ms.R, which must be run first
    aligned <- read_fa(file.path(out_dir, "prank_seeds", paste0("selected_", job$set, "_", job$lek, ".fa")))
  } else {
    input <- file.path(seq_dir, paste0("song_seqs_", job$set, "_", job$lek, ".fa"))
    matrix_file <- file.path(out_dir, paste0(job$strategy, "_cost_matrix.txt"))
    tmp <- file.path(work, basename(output))
    log <- file.path(out_dir, paste0(job$set, "_", job$lek, "_", job$strategy, "_mafft_output.txt"))
    system(paste(mafft, "--localpair --maxiterate 50000 --textmatrix", matrix_file, input, ">", tmp, "2>", log))
    aligned <- read_fa(tmp)
  }
  write_fa(aligned, output)
  output
}

# remaining years: each one takes the aligned row of its song type in the early alignment
song_type <- function(taxa) sub("-[0-9]{4}$", "", taxa)

add_years <- function(job) {
  early <- read_fa(file.path(out_dir, paste0(job$lek, "_", job$strategy, "_early_", job$fossils, "_alignment.fa")))
  output <- file.path(out_dir, paste0(job$lek, "_", job$strategy, "_all_", job$fossils, "_alignment.fa"))
  suffix <- if (job$strategy == "prank") "_prank.fa" else ".fa"
  rest <- read_fa(file.path(seq_dir, paste0("song_seqs_rest_", job$fossils, "_", job$lek, suffix)))
  rows <- setNames(early, song_type(names(early)))
  added <- setNames(rows[song_type(names(rest))], names(rest))
  if (anyNA(added)) stop("song type missing from the early alignment: ", job$lek, " ", job$fossils)
  all_rows <- c(early, added)
  write_fa(all_rows[order(match(song_type(names(all_rows)), song_type(names(early))), names(all_rows))], output)
  output
}

# first and last year of each song type
bound_taxa <- function(taxa) {
  year <- as.numeric(sub("^.*-", "", taxa))
  type <- song_type(taxa)
  taxa[year == ave(year, type, FUN = min) | year == ave(year, type, FUN = max)]
}

write_bounds <- function(job) {
  all_rows <- read_fa(file.path(out_dir, paste0(job$lek, "_", job$strategy, "_all_", job$fossils, "_alignment.fa")))
  output <- file.path(out_dir, paste0(job$lek, "_", job$strategy, "_bounds_", job$fossils, "_alignment.fa"))
  write_fa(all_rows[bound_taxa(names(all_rows))], output)
  output
}

fa <- list.files(seq_dir, pattern = "^song_seqs_(early_new|early_old|no\\.fossils)_[A-Z0-9]+\\.fa$")
m <- regmatches(fa, regexec("^song_seqs_(.+)_([A-Z0-9]+)\\.fa$", fa))
direct <- do.call(rbind, lapply(c("optimal", "all.equal", "prank"), function(s)
  data.frame(set = sapply(m, `[`, 2), lek = sapply(m, `[`, 3), strategy = s, stringsAsFactors = FALSE)))
early <- direct[grepl("^early", direct$set), ]
early$fossils <- sub("early_", "", early$set)

t0 <- Sys.time()
run <- function(jobs, f) {
  o <- parallel::mclapply(seq_len(nrow(jobs)), mc.cores = cores, mc.preschedule = FALSE, function(i) f(jobs[i, ]))
  bad <- vapply(o, inherits, TRUE, "try-error")
  if (any(bad)) stop("failed jobs: ", paste(which(bad), collapse = ", "), "\n", paste(o[bad], collapse = "\n"))
}
run(direct, align_direct)
run(early, add_years)
run(early, write_bounds)

# check every alignment against its input
check <- function(file) {
  p <- regmatches(basename(file), regexec("^([A-Z0-9]+)_(optimal|all\\.equal|prank)_(.+)_alignment\\.fa$", basename(file)))[[1]]
  lek <- p[2]; strategy <- p[3]; set <- p[4]
  suffix <- if (strategy == "prank") "_prank.fa" else ".fa"
  inputs <- if (grepl("^(all|bounds)_", set)) {
    fossils <- sub("^(all|bounds)_", "", set)
    c(read_fa(file.path(seq_dir, paste0("song_seqs_early_", fossils, "_", lek, suffix))),
      read_fa(file.path(seq_dir, paste0("song_seqs_rest_", fossils, "_", lek, suffix))))
  } else read_fa(file.path(seq_dir, paste0("song_seqs_", set, "_", lek, suffix)))
  if (grepl("^bounds_", set)) inputs <- inputs[bound_taxa(names(inputs))]

  aln <- read_fa(file)
  stripped <- gsub("-", "", aln)
  same_taxa <- setequal(names(aln), names(inputs))
  same_seqs <- same_taxa && all(stripped[names(inputs)] == inputs)
  aln_len <- nchar(aln[1])

  # all years of a song type have the same sequence, so they should have the same aligned row
  type <- sub("-[0-9]{4}$", "", names(aln))
  rows_per_type <- tapply(aln, type, function(z) length(unique(z)))

  log <- if (strategy == "prank") character() else {
    l <- file.path(out_dir, paste0(set, "_", lek, "_", strategy, "_mafft_output.txt"))
    if (file.exists(l)) readLines(l) else character()
  }
  data.frame(lek, strategy, set, taxa = length(aln), aln_len, same_length = all(nchar(aln) == aln_len),
             gaps = aln_len - max(nchar(inputs)), gap_prop = sum(nchar(aln) - nchar(stripped)) / sum(nchar(aln)),
             converged = if (strategy == "prank" || grepl("^(all|bounds)_", set)) NA else any(grepl("Converged", log)),
             same_taxa, same_seqs, types_with_different_rows = sum(rows_per_type > 1),
             stringsAsFactors = FALSE)
}

files <- list.files(out_dir, pattern = "_alignment\\.fa$", full.names = TRUE)
res <- do.call(rbind, lapply(files, check))
write.csv(res, file.path(out_dir, "alignment_results.csv"), row.names = FALSE)

unlink(work, recursive = TRUE)
cat("alignments:", nrow(res), "\n")
print(table(res$strategy, res$set))
cat("wrong taxa:", sum(!res$same_taxa), " wrong sequences:", sum(!res$same_seqs), " uneven rows:", sum(!res$same_length),
    " alignments with years of a song type in different rows:", sum(res$types_with_different_rows > 0), "\n")
print(Sys.time() - t0)
