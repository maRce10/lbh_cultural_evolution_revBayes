# Marginal likelihoods of the BEAST2 tree priors for the 5 ms pilot (bounds, relaxed clock):
# SA (from make_ml.R), sRanges (SR) and sRanges with single-occurrence ranges (SRsingle, the calibration).
# Same estimator as summarise_ml.R (stepping stone and path sampling from the likelihood log of each power).
#
# Run from the project root, after run_ml_sranges.sh:
#   Rscript scripts/revbayes_beast2_comparison/ml/summarise_ml_sranges.R

out <- "./output/revbayes_beast2_comparison/ml"
work_b2 <- file.path(Sys.getenv("HOME"), "lbh_ml_work")
n_powers <- 100
alpha <- 0.2
b2_burnin <- 0.6

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

rows <- list()
for (l in c("HC1", "BR1")) for (prior in c("SA", "SR", "SRsingle")) for (r in 1:3) {
  run <- if (prior == "SA") sprintf("%s_bounds_relaxed_r%d", l, r) else sprintf("%s_bounds_%s_relaxed_r%d", l, prior, r)
  screen <- file.path(out, "beast2", paste0(run, "_screen.txt"))
  if (!file.exists(screen) || !any(grepl("^marginal L estimate = ", readLines(screen)))) next
  b <- read_b2(file.path(work_b2, run))
  if (is.null(b)) next
  rows[[length(rows) + 1]] <- data.frame(lek = l, tree_prior = prior, replicate = r, t(ss_ps(b$p, b$L)))
}
ml <- do.call(rbind, rows)
write.csv(ml, file.path(out, "ml_sranges.csv"), row.names = FALSE)
print(ml, digits = 5)
print(aggregate(cbind(stepping_stone, path_sampling) ~ lek + tree_prior, ml, function(x) round(c(mean = mean(x), sd = sd(x)), 2)))
