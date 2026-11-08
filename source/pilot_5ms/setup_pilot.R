# Set up the RevBayes pilot of the 5 ms sequences: HC1, BR1 and SUR, MAFFT-optimal, new records, all and bounds,
# relaxed clock.
#
# Run from the project root:
#   Rscript source/pilot_5ms/setup_pilot.R [generations] [burnin]

args <- commandArgs(TRUE)
generations <- if (length(args) > 0) as.integer(args[1]) else NA   # NA: same as the first version of the study
burnin <- if (length(args) > 1) as.integer(args[2]) else NA

leks <- c("HC1", "BR1", "SUR")
datasets <- c("all", "bounds")
out <- "./data/processed/pilot_5ms"
rev_dir <- "./source/pilot_5ms"
for (d in file.path(out, c("nexus", "taxa"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

read_fa <- function(file) {
  a <- readLines(file)
  starts <- grep("^>", a)
  s <- mapply(function(i, j) paste(a[(i + 1):j], collapse = ""), starts, c(starts[-1] - 1, length(a)))
  setNames(as.vector(s), sub("^>", "", a[starts]))
}
write_fa <- function(s, file) writeLines(as.vector(rbind(paste0(">", names(s)), s)), file)
song_type <- function(x) sub("-[0-9]{4}$", "", x)
taxon_year <- function(x) as.numeric(sub("^.*-", "", x))

write_nexus <- function(aln, file) {
  writeLines(c("#NEXUS", "Begin data;", sprintf("dimensions ntax=%d nchar=%d;", length(aln), nchar(aln[1])),
               'format datatype=standard gap=- missing=? symbols="dfpsu";', "matrix", "",
               paste0(names(aln), "\t", aln), ";", "end;"), file)
}

summary_rows <- list()
for (l in leks) {
  all_5ms <- read_fa(file.path("./data/processed/mafft_alignments_5ms", paste0(l, "_optimal_all_new_alignment.fa")))
  bounds_5ms <- read_fa(file.path("./data/processed/mafft_alignments_5ms", paste0(l, "_optimal_bounds_new_alignment.fa")))

  sets <- list("5ms_all" = all_5ms, "5ms_bounds" = bounds_5ms)

  for (s in names(sets)) {
    aln <- sets[[s]]
    if (anyNA(aln)) stop("missing song type in ", l, " ", s)
    id <- paste(l, s, sep = "_")
    nexus <- file.path(out, "nexus", paste0(id, ".nex"))
    taxa <- file.path(out, "taxa", paste0(id, ".tsv"))
    write_nexus(aln, nexus)
    last <- max(taxon_year(names(aln)))
    age <- last - taxon_year(names(aln))
    writeLines(c("taxon\tage", paste(names(aln), age, sep = "\t")), taxa)

    extant <- names(aln)[age == 0]
    long_run <- l == "SUR"                       # as in the first version: longer runs for SUR
    gens <- if (is.na(generations)) if (long_run) 300000 else 150000 else generations
    burn <- if (is.na(burnin)) if (long_run) 100000 else 50000 else burnin
    tune <- max(10, burn %/% 250)
    writeLines(c(
      sprintf('alignment = "%s"', sub("^\\./", "./", nexus)),
      sprintf('fossils = "%s"', taxa),
      sprintf('base_str = "./output/pilot_5ms/%s"', id),
      sprintf('clade_extant = clade(%s)', paste0('"', extant, '"', collapse = ", ")),
      sprintf("burnin_generations = %d", burn),
      sprintf("run_generations = %d", gens),
      sprintf("tuning_interval = %d", tune),
      'source("./source/pilot_5ms/control.Rev")'), file.path(rev_dir, paste0("run_", id, ".Rev")))

    summary_rows[[id]] <- data.frame(lek = l, dataset = sub(".*_", "", s), taxa = length(aln),
                                     characters = nchar(aln[1]), gap_prop = sum(nchar(gsub("[^-]", "", aln))) / sum(nchar(aln)),
                                     taxa_last_year = length(extant), burnin = burn, generations = gens)
  }
}
info <- do.call(rbind, summary_rows)
rownames(info) <- NULL
write.csv(info, file.path(out, "pilot_datasets.csv"), row.names = FALSE)
print(info, digits = 2)
