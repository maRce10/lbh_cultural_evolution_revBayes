# Improvement experiments for 6-state classification of 20 equal-length segments
library(xgboost)
setwd("/home/m/Dropbox/Projects/lbh_cultural_evolution_revBayes")

X <- read.csv("./data/processed/segment_classifier_20/segment_features_20.csv")
contours <- read.csv("./data/raw/segment_classification/seltailor_output.csv")
ff_cols <- grep("^ff", names(contours))
X <- X[order(X$ID, X$segment), ]
classes6 <- c("D", "U", "P", "F", "M", "S")
y6 <- match(X$Category, classes6) - 1
yb <- as.integer(X$Category %in% c("F", "M", "S"))
non_feat <- c("ID", "segment", "start", "end", "Category", "purity", "lek.song.type", "year", "paper", "paper_used")
feat <- c("segment", setdiff(names(X), non_feat))
for (v in feat) X[[v]][is.na(X[[v]])] <- -999

params <- list(eta = 0.05, max_depth = 6, subsample = 0.8, colsample_bytree = 0.5, min_child_weight = 2)
cw6 <- as.numeric(sqrt(max(table(y6)) / table(y6)))[y6 + 1]
cwb <- ifelse(yb == 1, sqrt(sum(yb == 0) / sum(yb == 1)), 1)

fit <- function(rows, M, lab, nc, w, nr = 500) {
  obj <- if (nc == 2) list(objective = "binary:logistic") else list(objective = "multi:softprob", num_class = nc)
  xgb.train(params = c(params, obj, nthread = 16), data = xgb.DMatrix(as.matrix(M[rows, , drop = FALSE]), label = lab[rows], weight = w[rows]),
            nrounds = nr, verbose = 0)
}
pred <- function(m, M, rows, nc) {
  p <- predict(m, as.matrix(M[rows, , drop = FALSE]))
  if (nc == 2) p else matrix(p, ncol = nc, byrow = TRUE)
}

# contour features over a time window (same as in the qmd, reduced)
cf <- function(tt, ff, a, b) {
  n <- max(8, sum(tt >= a & tt <= b)); ti <- seq(a, b, length.out = n); fi <- approx(tt, ff, xout = ti, rule = 2)$y
  dur <- b - a; lin <- lm(fi ~ ti); d1 <- diff(fi); sig <- sign(d1[abs(d1) > 0.05])
  n_ext <- if (length(sig) > 1) sum(diff(sig) != 0) else 0
  c(run_dur = dur * 1000, run_osc = n_ext / 2 / dur, run_period = if (n_ext > 0) 1000 / (n_ext / 2 / dur) else 100,
    run_slope = unname(coef(lin)[2]), run_range = diff(range(fi)), run_resid = sd(residuals(lin)), run_fmean = mean(fi))
}
song_ct <- lapply(split(contours, contours$song.type), function(ct) list(tt = seq(ct$start, ct$end, length.out = length(ff_cols)), ff = unlist(ct[, ff_cols])))

# run-level features given a predicted trill indicator
run_features <- function(is_trill) {
  out <- matrix(NA, nrow(X), 9)
  for (id in unique(X$ID)) {
    idx <- which(X$ID == id); s <- song_ct[[id]]
    r <- rle(is_trill[idx]); e <- cumsum(r$lengths); b <- e - r$lengths + 1
    for (k in seq_along(r$values)) {
      ii <- idx[b[k]:e[k]]
      f <- cf(s$tt, s$ff, min(X$start[ii]), max(X$end[ii]))
      out[ii, ] <- matrix(c(f, NA, NA), nrow = length(ii), ncol = 9, byrow = TRUE)
      out[ii, 8] <- seq_along(ii) / length(ii)      # relative position within run
      out[ii, 9] <- length(ii)                       # run length (segments)
    }
  }
  colnames(out) <- c("run_dur", "run_osc", "run_period", "run_slope", "run_range", "run_resid", "run_fmean", "run_relpos", "run_len")
  out
}

# neighbour probabilities (+-1, +-2) within songs
neigh <- function(P) {
  n <- nrow(P); out <- NULL
  for (k in c(-2, -1, 1, 2)) {
    idx <- seq_len(n) + k; ok <- idx >= 1 & idx <= n
    ok[ok] <- X$ID[idx[ok]] == X$ID[ok]
    Q <- matrix(-1, n, ncol(P)); Q[ok, ] <- P[idx[ok], ]
    colnames(Q) <- paste0("nb", k, "_", seq_len(ncol(P))); out <- cbind(out, Q)
  }
  out
}

# Viterbi decoding per song using transition matrix from training songs
viterbi_decode <- function(P, train_rows, test_rows) {
  tr_lab <- y6[train_rows]; tr_id <- X$ID[train_rows]
  Tm <- matrix(1, 6, 6)
  for (i in seq_along(tr_lab)[-1]) if (tr_id[i] == tr_id[i - 1]) Tm[tr_lab[i - 1] + 1, tr_lab[i] + 1] <- Tm[tr_lab[i - 1] + 1, tr_lab[i] + 1] + 1
  Tm <- log(Tm / rowSums(Tm))
  prior <- tabulate(tr_lab + 1, 6) / length(tr_lab)
  first <- tabulate(tr_lab[X$segment[train_rows] == 1] + 1, 6) + 1; first <- log(first / sum(first))
  out <- rep(NA, nrow(X))
  for (id in unique(X$ID[test_rows])) {
    ii <- test_rows[X$ID[test_rows] == id]; ii <- ii[order(X$segment[ii])]
    E <- log(pmax(P[ii, , drop = FALSE], 1e-6)) - matrix(log(prior), length(ii), 6, byrow = TRUE) * 0.5
    V <- matrix(-Inf, length(ii), 6); B <- matrix(0, length(ii), 6)
    V[1, ] <- first + E[1, ]
    for (t in seq_along(ii)[-1]) for (j in 1:6) { s <- V[t - 1, ] + Tm[, j]; B[t, j] <- which.max(s); V[t, j] <- max(s) + E[t, j] }
    path <- integer(length(ii)); path[length(ii)] <- which.max(V[length(ii), ])
    for (t in rev(seq_along(ii)[-1])) path[t - 1] <- B[t, path[t]]
    out[ii] <- path - 1
  }
  out
}

acc <- function(p) c(acc = mean(p == y6), acc_pure = mean((p == y6)[X$purity >= 0.8]))

res <- list()
Xf <- X[, feat]
for (r in 1:2) {
  set.seed(r); types <- unique(X$lek.song.type)
  fold <- setNames(sample(rep(1:5, length.out = length(types))), types)[X$lek.song.type]
  P_base <- P_run <- P_stack <- P_pur <- matrix(NA, nrow(X), 6)
  vit_base <- vit_stack <- rep(NA, nrow(X))
  for (f in 1:5) {
    tr <- which(fold != f); te <- which(fold == f)
    # 1 baseline
    m <- fit(tr, Xf, y6, 6, cw6); P_base[te, ] <- pred(m, Xf, te, 6)
    # 3 purity weighting
    m <- fit(tr, Xf, y6, 6, cw6 * X$purity); P_pur[te, ] <- pred(m, Xf, te, 6)
    # inner CV for binary (run features) and 6-state probs (stacking) on training rows
    itypes <- unique(X$lek.song.type[tr])
    ifold <- setNames(sample(rep(1:4, length.out = length(itypes))), itypes)[X$lek.song.type[tr]]
    pb <- rep(NA, nrow(X)); P_in <- matrix(NA, nrow(X), 6)
    for (g in 1:4) {
      a <- tr[ifold != g]; b <- tr[ifold == g]
      pb[b] <- pred(fit(a, Xf, yb, 2, cwb), Xf, b, 2)
      P_in[b, ] <- pred(fit(a, Xf, y6, 6, cw6), Xf, b, 6)
    }
    pb[te] <- pred(fit(tr, Xf, yb, 2, cwb), Xf, te, 2)
    P_in[te, ] <- P_base[te, ]
    # 2 run-level features (runs defined by predicted tone/trill)
    RF <- run_features(pb > 0.5); Xr <- cbind(Xf, RF, p_trill = pb); Xr[is.na(Xr)] <- -999
    m <- fit(tr, Xr, y6, 6, cw6); P_run[te, ] <- pred(m, Xr, te, 6)
    # 4 stacking: run features + own and neighbour probabilities
    Xs <- cbind(Xr, own = P_in, neigh(P_in)); Xs[is.na(Xs)] <- -999
    m <- fit(tr, Xs, y6, 6, cw6); P_stack[te, ] <- pred(m, Xs, te, 6)
    # 5 Viterbi decoding
    vit_base[te] <- viterbi_decode(P_base, tr, te)[te]
    vit_stack[te] <- viterbi_decode(P_stack, tr, te)[te]
  }
  res[[r]] <- rbind(
    "1 baseline" = acc(max.col(P_base) - 1),
    "2 + run-level features" = acc(max.col(P_run) - 1),
    "3 baseline + purity weights" = acc(max.col(P_pur) - 1),
    "4 + run features + neighbour stacking" = acc(max.col(P_stack) - 1),
    "5a baseline + Viterbi" = acc(vit_base),
    "5b stacking + Viterbi" = acc(vit_stack))
  cat("repetition", r, "\n"); print(round(res[[r]], 3))
  if (r == 1) saveRDS(list(P_base = P_base, P_stack = P_stack, vit_stack = vit_stack, fold = fold), "/tmp/claude-1000/-home-m-Dropbox/c31c64a2-d619-45c1-a333-94f805bdb4ff/scratchpad/improve_20seg_oof.RDS")
}
cat("\nmean over repetitions:\n"); print(round(Reduce(`+`, res) / length(res), 3))
