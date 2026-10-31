# Train and evaluate segment category classifiers on contour features
# Cross-validation holds out whole lek song types (all years of e.g. BR1-A1 together)

library(randomForest)
library(xgboost)

set.seed(123)
X <- read.csv("./data/processed/segment_classifier/segment_features.csv")
X$Category <- factor(X$Category)
id_cols <- c("ID", "lek.song.type.year", "lek.song.type", "Category", "start", "end", "song_f_mean")
feat <- setdiff(names(X), id_cols)
for (v in feat) X[[v]][is.na(X[[v]])] <- -999   # first/last segments have no neighbour

# folds by song type
k <- 5
types <- unique(X$lek.song.type)
fold_of_type <- setNames(sample(rep(seq_len(k), length.out = length(types))), types)
X$fold <- fold_of_type[X$lek.song.type]

classes <- levels(X$Category)
p_rf <- p_xgb <- matrix(NA, nrow(X), length(classes), dimnames = list(NULL, classes))
for (f in seq_len(k)) {
  tr <- X$fold != f; te <- X$fold == f
  rf <- randomForest(x = X[tr, feat], y = X$Category[tr], ntree = 1000)
  p_rf[te, ] <- predict(rf, X[te, feat], type = "prob")
  xg <- xgboost(data = as.matrix(X[tr, feat]), label = as.integer(X$Category[tr]) - 1,
                objective = "multi:softprob", num_class = length(classes), nrounds = 400,
                eta = 0.05, max_depth = 4, subsample = 0.8, colsample_bytree = 0.6,
                min_child_weight = 2, verbose = 0, nthread = 8)
  p_xgb[te, ] <- matrix(predict(xg, as.matrix(X[te, feat])), ncol = length(classes), byrow = TRUE)
}

evaluate <- function(p, name) {
  pred <- factor(classes[max.col(p)], levels = classes)
  cm <- table(true = X$Category, pred = pred)
  acc <- sum(diag(cm)) / sum(cm)
  pe <- sum(rowSums(cm) * colSums(cm)) / sum(cm)^2
  cat(sprintf("\n%s: accuracy = %.3f, kappa = %.3f\n", name, acc, (acc - pe) / (1 - pe)))
  print(cm)
  cat("per-class recall:\n"); print(round(diag(cm) / rowSums(cm), 2))
  invisible(pred)
}
pred_rf <- evaluate(p_rf, "Random forest")
pred_xgb <- evaluate(p_xgb, "xgboost")
pred_ens <- evaluate((p_rf + p_xgb) / 2, "Ensemble (RF + xgboost)")

# fit final model on all data and save
rf_all <- randomForest(x = X[, feat], y = X$Category, ntree = 1000, importance = TRUE)
saveRDS(list(model = rf_all, features = feat), "./data/processed/segment_classifier/rf_segment_classifier.RDS")
imp <- importance(rf_all, type = 1)
cat("\nTop features (mean decrease accuracy):\n"); print(round(head(imp[order(-imp[, 1]), , drop = FALSE], 15), 1))

out <- cbind(X[, c("ID", "lek.song.type", "start", "end", "Category")], pred = pred_ens,
             conf = apply((p_rf + p_xgb) / 2, 1, max))
write.csv(out, "./data/processed/segment_classifier/cv_predictions.csv", row.names = FALSE)
