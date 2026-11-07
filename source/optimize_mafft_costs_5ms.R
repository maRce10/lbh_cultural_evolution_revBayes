# Optimize the MAFFT scoring matrix for the 5 ms, 5-state song sequences.
#
# Follows source/Optimizing_alignment.Rmd: real songs should need fewer gaps than
# simulated songs built from the same runs of states, so the best scores are those that
# maximize the difference in gaps (simulated - real). Scores always increase:
# different category (trill vs. tone) < same category < same state.
#
# Run from the project root:
#   Rscript source/optimize_mafft_costs_5ms.R <n simulated sets> <cores> [n cost combinations, for a test]

args <- commandArgs(TRUE)
n_sets <- if (length(args) > 0) as.integer(args[1]) else 100
cores <- if (length(args) > 1) as.integer(args[2]) else 28
n_test <- if (length(args) > 2) as.integer(args[3]) else NA

# mafft v7.526 from bioconda: mamba create -n mafft -c conda-forge -c bioconda mafft
mafft <- path.expand("~/mambaforge/envs/mafft/bin/mafft")
if (!file.exists(mafft)) mafft <- Sys.which("mafft")
if (!nzchar(mafft)) stop("mafft not found")
out_dir <- "./data/processed/mafft_optimization_5ms"
res_dir <- file.path(out_dir, "by_cost")
work <- tempfile("mafft_opt_")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work)

# real songs: leks with more than 5 song types
dat <- read.csv("./data/processed/song_5ms_sequences_by_song_type_all_leks.csv", stringsAsFactors = FALSE)
tb <- table(dat$lek)
leks <- names(tb)[tb > 5]
dat <- dat[dat$lek %in% leks, ]

# lexicon: runs of the same state in the real songs
runs <- function(x) { r <- rle(strsplit(x, "")[[1]]); strrep(r$values, r$lengths) }
lexic <- unlist(lapply(dat$sequence, runs), use.names = FALSE)

# a simulated song with the same length as a real song, with no two consecutive runs of the same state
sim_song <- function(len) {
  repeat {
    song <- ""
    while (nchar(song) < len) {
      left <- len - nchar(song)
      last <- substr(song, nchar(song), nchar(song))
      cand <- lexic[nchar(lexic) <= left & substr(lexic, 1, 1) != last]
      if (!length(cand)) break
      song <- paste0(song, sample(cand, 1))
    }
    if (nchar(song) == len) return(song)
  }
}

# each simulated lek has the same number of songs and the same song lengths as the real lek
set.seed(10)
sim_sets <- lapply(seq_len(n_sets), function(s)
  lapply(setNames(leks, leks), function(l) vapply(nchar(dat$sequence[dat$lek == l]), sim_song, "")))
saveRDS(sim_sets, file.path(out_dir, "simulated_song_sets.RDS"))
writeLines(lexic, file.path(out_dir, "song_state_lexicon.txt"))

write_fa <- function(seqs, file) writeLines(as.vector(rbind(paste0(">s", seq_along(seqs)), seqs)), file)

inputs <- list()
for (l in leks) {
  f <- file.path(work, paste0("real_", l, ".fa")); write_fa(dat$sequence[dat$lek == l], f)
  inputs[[f]] <- data.frame(file = f, lek = l, data.type = "real", set = 0, max_len = max(nchar(dat$sequence[dat$lek == l])))
  for (s in seq_len(n_sets)) {
    f <- file.path(work, paste0("sim", s, "_", l, ".fa")); write_fa(sim_sets[[s]][[l]], f)
    inputs[[f]] <- data.frame(file = f, lek = l, data.type = "simulated", set = s, max_len = max(nchar(sim_sets[[s]][[l]])))
  }
}
inputs <- do.call(rbind, inputs)

# scores: always increasing, with a positive score for the same state
costs <- seq(-2, 2, by = 0.25)
grid <- as.data.frame(t(combn(costs, 3)))
names(grid) <- c("different.category", "same.category", "same.state")
grid <- grid[grid$same.state > 0, ]
grid$cost.vals <- paste(grid$different.category, grid$same.category, grid$same.state, sep = "_")
if (!is.na(n_test)) grid <- grid[round(seq(1, nrow(grid), length.out = n_test)), ]

# each state is written with its own ascii code
states <- c("d", "f", "p", "s", "u")
hex <- c(d = "0x64", f = "0x66", p = "0x70", s = "0x73", u = "0x75")
trills <- c("f", "s")
pairs <- expand.grid(a = states, b = states, stringsAsFactors = FALSE)

write_matrix <- function(g, file) {
  v <- ifelse(pairs$a == pairs$b, g$same.state,
              ifelse((pairs$a %in% trills) == (pairs$b %in% trills), g$same.category, g$different.category))
  writeLines(paste(hex[pairs$a], hex[pairs$b], v), file)
}

align <- function(input, matrix_file) {
  tag <- tempfile(tmpdir = work)
  on.exit(unlink(paste0(tag, c(".fa", ".log"))))
  system(paste(mafft, "--localpair --maxiterate 10000 --textmatrix", matrix_file, input$file,
               ">", paste0(tag, ".fa"), "2>", paste0(tag, ".log")))
  a <- readLines(paste0(tag, ".fa"))
  if (!length(a)) return(data.frame(aln_len = NA, gaps = NA, gap_prop = NA, convg = "N"))
  starts <- grep("^>", a)
  aln <- mapply(function(i, j) paste(a[(i + 1):j], collapse = ""), starts, c(starts[-1] - 1, length(a)))
  n_gap <- sum(nchar(aln) - nchar(gsub("-", "", aln)))
  data.frame(aln_len = nchar(aln[1]), gaps = nchar(aln[1]) - input$max_len, gap_prop = n_gap / sum(nchar(aln)),
             convg = ifelse(any(grepl("Converged", readLines(paste0(tag, ".log")))), "Y", "N"))
}

t0 <- Sys.time()
done <- parallel::mclapply(seq_len(nrow(grid)), mc.cores = cores, mc.preschedule = FALSE, function(k) {
  g <- grid[k, ]
  out <- file.path(res_dir, paste0(g$cost.vals, ".csv"))
  if (file.exists(out)) return(TRUE)
  m <- file.path(work, paste0("matrix_", k, ".txt"))
  write_matrix(g, m)
  res <- do.call(rbind, lapply(seq_len(nrow(inputs)), function(i) cbind(inputs[i, c("lek", "data.type", "set")], align(inputs[i, ], m))))
  res$cost.vals <- g$cost.vals
  write.csv(res, paste0(out, ".tmp"), row.names = FALSE)
  file.rename(paste0(out, ".tmp"), out)
  unlink(m)
  TRUE
})
unlink(work, recursive = TRUE)
cat("cost combinations:", nrow(grid), " alignments per combination:", nrow(inputs),
    " done:", sum(unlist(lapply(done, isTRUE))), "\n")
print(Sys.time() - t0)
