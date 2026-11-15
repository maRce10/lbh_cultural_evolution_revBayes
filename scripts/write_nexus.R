# Convert the 5 ms song alignments to nexus files for RevBayes, with the format used in the previous study,
# and write the matching taxon files (taxon and age in years before the last year of sampling of the lek).
# Only the sets to be analyzed are written (early is only the base alignment):
#   dataset: all (every year of each song type) or bounds (first and last year of each song type)
#   historical records: full (including recordings before 2000), recent (after 2000) or none (only the song
#   types of the last year of sampling)
#
# File names:
#   nexus: data/processed/nexus/<lek>_<strategy>_<dataset>_<historical>_alignment.nex, or
#          <lek>_<strategy>_none_alignment.nex
#   taxa:  data/processed/taxa/<lek>_<dataset>_<historical>_taxa.tsv, or <lek>_none_taxa.tsv
#          (the same for the three alignment strategies)
#
# Run from the project root, after scripts/align_song_types.R:
#   Rscript scripts/write_nexus.R

in_dir <- "./data/processed/alignments"
out_dir <- "./data/processed/nexus"
taxa_dir <- "./data/processed/taxa"
dir.create(out_dir, showWarnings = FALSE)
dir.create(taxa_dir, showWarnings = FALSE)

read_fa <- function(file) {
  a <- readLines(file)
  starts <- grep("^>", a)
  s <- mapply(function(i, j) paste(a[(i + 1):j], collapse = ""), starts, c(starts[-1] - 1, length(a)))
  setNames(as.vector(s), sub("^>", "", a[starts]))
}

files <- list.files(in_dir, pattern = "_(all_recent|all_full|bounds_recent|bounds_full|none)_alignment\\.fa$", full.names = TRUE)

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
  p <- regmatches(basename(f), regexec("^([A-Z0-9]+)_(optimal|agnostic|prank)_(.+)_alignment\\.fa$", basename(f)))[[1]]
  data.frame(file = basename(f), lek = p[2], strategy = p[3], set = p[4], ntax = length(aln), nchar = nchar(aln[1]),
             taxa = paste(sort(names(aln)), collapse = " "))
}))

# taxon files: one per lek and set, checking that the three strategies have the same taxa
year <- function(taxa) as.numeric(sub("^.*-", "", taxa))
last_year <- tapply(unlist(lapply(strsplit(res$taxa, " "), year)), rep(res$lek, lengths(strsplit(res$taxa, " "))), max)
for (k in split(res, list(res$lek, res$set), drop = TRUE)) {
  if (length(unique(k$taxa)) != 1) stop("different taxa among strategies: ", k$lek[1], " ", k$set[1])
  taxa <- strsplit(k$taxa[1], " ")[[1]]
  writeLines(c("taxon\tage", paste(taxa, last_year[[k$lek[1]]] - year(taxa), sep = "\t")),
             file.path(taxa_dir, paste0(k$lek[1], "_", k$set[1], "_taxa.tsv")))
}

cat("nexus files:", nrow(res), " taxon files:", length(list.files(taxa_dir, pattern = "_taxa\\.tsv$")), "\n")
print(table(res$strategy, res$set))
