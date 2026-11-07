# Convert the 5 ms song alignments to nexus files for RevBayes, with the format used in the previous study.
# Only the sets to be analyzed are written: all, bounds and no.fossils (early is only the base alignment).
#
# Run from the project root, after source/align_song_types_5ms.R:
#   Rscript source/write_nexus_5ms.R

in_dir <- "./data/processed/mafft_alignments_5ms"
out_dir <- "./data/processed/nexus_5ms"
dir.create(out_dir, showWarnings = FALSE)

read_fa <- function(file) {
  a <- readLines(file)
  starts <- grep("^>", a)
  s <- mapply(function(i, j) paste(a[(i + 1):j], collapse = ""), starts, c(starts[-1] - 1, length(a)))
  setNames(as.vector(s), sub("^>", "", a[starts]))
}

files <- list.files(in_dir, pattern = "_(all_new|all_old|bounds_new|bounds_old|no\\.fossils)_alignment\\.fa$", full.names = TRUE)

res <- do.call(rbind, lapply(files, function(f) {
  aln <- read_fa(f)
  # PRANK alignments use nucleotide letters: back to song states
  if (grepl("_prank_", basename(f))) aln <- setNames(chartr("CTAGR", "fsudp", aln), names(aln))

  if (any(grepl("[^dfpsu-]", aln))) stop("unexpected symbols in ", f)
  if (length(unique(nchar(aln))) != 1) stop("rows of different length in ", f)
  if (anyDuplicated(names(aln))) stop("repeated taxa in ", f)

  writeLines(c("#NEXUS", "Begin data;", sprintf("dimensions ntax=%d nchar=%d;", length(aln), nchar(aln[1])),
               'format datatype=standard gap=- missing=? symbols="dfpsu";', "matrix", "",
               paste0(names(aln), "\t", aln), ";", "end;"),
             file.path(out_dir, sub("\\.fa$", ".nex", basename(f))))
  data.frame(file = basename(f), ntax = length(aln), nchar = nchar(aln[1]))
}))

cat("nexus files:", nrow(res), "\n")
print(table(sub("^[A-Z0-9]+_([a-z.]+)_.*$", "\\1", res$file), sub("^[A-Z0-9]+_[a-z.]+_(.*)_alignment.*$", "\\1", res$file)))
