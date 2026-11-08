# Collect the marginal likelihood estimates of the 5 ms pilot (RevBayes and BEAST2) into one table.
#
# Stepping stone (Xie et al. 2011) and path sampling (trapezoid rule) are computed here with the same code for
# both programs, from the log-likelihood samples at each power: the RevBayes .out files and the likelihood.log
# of each BEAST2 path sampler step. For RevBayes this reproduces RevBayes' own estimates.
# The estimate printed by the BEAST2 path sampler (PathSampleAnalyser) is kept for reference as beast2_reported:
# it does not match its own samples (it is off by tens of log units), so it is not used.
#
# Run from the project root:
#   Rscript source/pilot_5ms/ml/summarise_ml.R

out <- "./output/pilot_5ms/ml"
work_b2 <- file.path(Sys.getenv("HOME"), "lbh_ml_work")
n_powers <- 100
alpha <- 0.2
b2_burnin <- 0.6   # same as burnInPercentage in the BEAST2 XML

ss_ps <- function(powers, L) {
  o <- order(powers); p <- powers[o]; L <- L[o]
  ss <- sum(sapply(2:length(p), function(k) {
    d <- p[k] - p[k - 1]; x <- L[[k - 1]]; m <- max(x)
    d * m + log(mean(exp(d * (x - m))))
  }))
  mu <- sapply(L, mean)
  ps <- sum(diff(p) * (mu[-1] + mu[-length(mu)]) / 2)
  c(stepping_stone = ss, path_sampling = ps)
}

read_rb <- function(f) {
  x <- read.table(f, header = TRUE)
  s <- split(x$likelihood, x$power)
  list(p = as.numeric(names(s)), L = unname(s))
}

read_b2 <- function(dir) {
  files <- file.path(dir, paste0("step", 0:(n_powers - 1)), "likelihood.log")
  if (!all(file.exists(files))) return(NULL)
  L <- lapply(files, function(f) {
    x <- read.table(f, header = TRUE, comment.char = "#")
    n <- nrow(x)
    x[ceiling(b2_burnin * n):n, 2]
  })
  list(p = (((n_powers - 1):0) / (n_powers - 1))^(1 / alpha), L = L)
}

runs <- sub("\\.Rev$", "", list.files("./source/pilot_5ms/ml/revbayes", pattern = "\\.Rev$"))
rows <- list()
for (run in runs) {
  f <- file.path(out, "revbayes", paste0(run, ".out"))
  if (file.exists(f) && file.exists(file.path(out, "revbayes", paste0(run, "_ml.txt")))) {
    v <- read_rb(f)
    rows[[length(rows) + 1]] <- data.frame(program = "RevBayes", run = run, t(ss_ps(v$p, v$L)), beast2_reported = NA)
  }
  screen <- file.path(out, "beast2", paste0(run, "_screen.txt"))
  rep_line <- if (file.exists(screen)) grep("^marginal L estimate = ", readLines(screen), value = TRUE) else character()
  if (length(rep_line) > 0) {
    b <- read_b2(file.path(work_b2, run))
    if (!is.null(b))
      rows[[length(rows) + 1]] <- data.frame(program = "BEAST2", run = run, t(ss_ps(b$p, b$L)),
                                             beast2_reported = as.numeric(sub("^marginal L estimate = ", "", tail(rep_line, 1))))
  }
}

ml <- do.call(rbind, rows)
parts <- do.call(rbind, strsplit(ml$run, "_"))
ml$lek <- parts[, 1]; ml$dataset <- parts[, 2]; ml$clock <- parts[, 3]; ml$replicate <- parts[, 4]
ml <- ml[order(ml$lek, ml$dataset, ml$clock, ml$program, ml$replicate),
         c("program", "lek", "dataset", "clock", "replicate", "stepping_stone", "path_sampling", "beast2_reported")]
rownames(ml) <- NULL
write.csv(ml, file.path(out, "ml_estimates.csv"), row.names = FALSE)
print(ml, digits = 5)
