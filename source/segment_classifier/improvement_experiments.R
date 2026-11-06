# Compare strategies to improve segment classification
# All variants use the same grouped folds (whole lek song types held out), repeated 3 times
# Output: data/processed/segment_classifier/improvement_experiments.csv

library(xgboost)

X <- read.csv("./data/processed/segment_classifier/segment_features_audio.csv")
X$Category <- factor(X$Category)
classes <- levels(X$Category); nc <- length(classes)
id_cols <- c("ID", "lek.song.type.year", "lek.song.type", "Category", "start", "end", "song_f_mean")
audio_feat <- c("af_period_ms", "af_n_peaks", "af_ac_period_ms", "af_ac_peak", "af_range",
                "env_ac_period_ms", "env_ac_peak", "env_depth")
contour_feat <- setdiff(names(X), c(id_cols, audio_feat))
for (v in c(contour_feat, audio_feat)) X[[v]][is.na(X[[v]])] <- -999
y <- as.integer(X$Category) - 1

params_default <- list(eta = 0.05, max_depth = 4, subsample = 0.8, colsample_bytree = 0.6, min_child_weight = 2)

fit_predict <- function(tr, te, feat, params = params_default, nrounds = 400, weights = NULL, Xm = X) {
  d <- xgb.DMatrix(as.matrix(Xm[tr, feat]), label = y[tr], weight = weights[tr])
  m <- xgb.train(params = c(params, objective = "multi:softprob", num_class = nc, nthread = 8),
                 data = d, nrounds = nrounds, verbose = 0)
  matrix(predict(m, as.matrix(Xm[te, feat])), ncol = nc, byrow = TRUE)
}

# probabilities of the previous / next segment as extra features
neighbour_probs <- function(P, rows) {
  out <- matrix(-1, length(rows), 2 * nc)
  for (k in seq_along(rows)) {
    i <- rows[k]
    if (i > 1 && X$ID[i - 1] == X$ID[i] && !is.na(P[i - 1, 1])) out[k, 1:nc] <- P[i - 1, ]
    if (i < nrow(X) && X$ID[i + 1] == X$ID[i] && !is.na(P[i + 1, 1])) out[k, nc + 1:nc] <- P[i + 1, ]
  }
  colnames(out) <- c(paste0("prev_p_", classes), paste0("next_p_", classes)); out
}

metrics <- function(p) {
  pred <- max.col(p) - 1
  cm <- table(factor(y, 0:(nc - 1)), factor(pred, 0:(nc - 1)))
  acc <- sum(diag(cm)) / sum(cm); pe <- sum(rowSums(cm) * colSums(cm)) / sum(cm)^2
  rec <- diag(cm) / rowSums(cm); prec <- diag(cm) / pmax(colSums(cm), 1)
  f1 <- ifelse(rec + prec > 0, 2 * rec * prec / (rec + prec), 0)
  c(accuracy = acc, kappa = (acc - pe) / (1 - pe), macro_F1 = mean(f1), recall_F = rec[[which(classes == "F")]],
    recall_M = rec[[which(classes == "M")]], recall_P = rec[[which(classes == "P")]])
}

class_w <- as.numeric(sqrt(max(table(y)) / table(y)))[y + 1]   # sqrt-balanced class weights
all_feat <- c(contour_feat, audio_feat)
variants <- list(
  "1 contour features (baseline)" = list(feat = contour_feat),
  "2 + audio trill features" = list(feat = all_feat),
  "3 + class weights" = list(feat = all_feat, weights = class_w),
  "4 + deeper trees (depth 6, 800 rounds, eta 0.03)" = list(feat = all_feat, weights = class_w,
       params = modifyList(params_default, list(max_depth = 6, eta = 0.03)), nrounds = 800)
)

res <- list()
for (r in 1:3) {
  set.seed(r)
  types <- unique(X$lek.song.type)
  fold <- setNames(sample(rep(1:5, length.out = length(types))), types)[X$lek.song.type]
  P <- lapply(variants, function(v) matrix(NA, nrow(X), nc))
  P_stack <- matrix(NA, nrow(X), nc)
  for (f in 1:5) {
    tr <- which(fold != f); te <- which(fold == f)
    for (v in names(variants)) {
      a <- variants[[v]]
      P[[v]][te, ] <- fit_predict(tr, te, a$feat, if (is.null(a$params)) params_default else a$params,
                                  if (is.null(a$nrounds)) 400 else a$nrounds, a$weights)
    }
    # 5 stacking: stage 1 = variant 3; stage-1 probabilities for the training rows come from an
    # inner grouped CV, so stage 2 never sees probabilities from models trained on its own rows
    inner_P <- matrix(NA, nrow(X), nc)
    inner_fold <- setNames(sample(rep(1:4, length.out = length(unique(X$lek.song.type[tr])))),
                           unique(X$lek.song.type[tr]))[X$lek.song.type[tr]]
    for (g in 1:4) inner_P[tr[inner_fold == g], ] <- fit_predict(tr[inner_fold != g], tr[inner_fold == g], all_feat, weights = class_w)
    inner_P[te, ] <- P[["3 + class weights"]][te, ]
    Xs <- cbind(X[, all_feat], neighbour_probs(inner_P, seq_len(nrow(X))), own = inner_P)
    P_stack[te, ] <- fit_predict(tr, te, names(Xs), weights = class_w, Xm = Xs)
  }
  P[["5 + stacking with neighbour predictions"]] <- P_stack
  res[[r]] <- do.call(rbind, lapply(names(P), function(v) data.frame(variant = v, repeat_ = r, t(metrics(P[[v]])))))
  cat("repeat", r, "done\n")
}
res <- do.call(rbind, res)
summ <- aggregate(. ~ variant, res[, -2], function(z) sprintf("%.3f ± %.3f", mean(z), sd(z)))
print(summ, right = FALSE)
write.csv(res, "./data/processed/segment_classifier/improvement_experiments.csv", row.names = FALSE)
