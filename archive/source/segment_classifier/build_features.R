# Build per-segment features from tailored frequency contours
# Inputs: data/raw/segment_classification/seltailor_output.csv (150-point ff contours per song, kHz)
#         data/raw/segment_classification/FrameworkContourDegradationv1.7.txt (manually delimited and labelled segments)
# Output: data/processed/segment_classifier/segment_features.csv

library(Rraven)

contours <- read.csv("./data/raw/segment_classification/seltailor_output.csv")
ff_cols <- grep("^ff", names(contours))

segs <- imp_raven(path = "./data/raw/segment_classification", files = "FrameworkContourDegradationv1.7.txt",
                  all.data = TRUE, warbler.format = TRUE)
segs <- segs[!is.na(segs$Category) & segs$Category != "", ]
segs <- segs[order(segs$ID, segs$start), ]

# features describing the contour shape within a time window
contour_features <- function(tt, ff, a, b, n_resample = 20) {
  inside <- tt >= a & tt <= b
  # interpolate so very short segments still get a usable number of points
  n <- max(8, sum(inside))
  ti <- seq(a, b, length.out = n)
  fi <- approx(tt, ff, xout = ti, rule = 2)$y
  dur <- b - a
  lin <- lm(fi ~ ti)
  res <- residuals(lin)
  ts <- (ti - a) / dur   # time rescaled to 0-1 to keep the quadratic fit well conditioned
  quad <- coef(lm(fi ~ ts + I(ts^2)))[3]
  d1 <- diff(fi)
  # local extrema (direction changes ignoring tiny wiggles below 50 Hz)
  sig <- sign(d1[abs(d1) > 0.05])
  n_ext <- if (length(sig) > 1) sum(diff(sig) != 0) else 0
  # trill pace: median interval between successive contour peaks (and troughs)
  fd <- fi - fitted(lin)
  pk <- which(diff(sign(diff(fd))) == -2) + 1
  tr <- which(diff(sign(diff(fd))) == 2) + 1
  pk <- pk[fd[pk] > 0.1]; tr <- tr[fd[tr] < -0.1]
  ivl <- c(diff(ti[pk]), diff(ti[tr]))
  period <- if (length(ivl)) median(ivl) else NA
  # autocorrelation-based period of the detrended contour
  ac <- if (n > 10 && sd(fd) > 0) acf(fd, lag.max = n - 2, plot = FALSE)$acf[-1] else NA
  ac_lag <- if (length(ac) > 2 && any(diff(sign(diff(ac))) == -2)) {
    lp <- which(diff(sign(diff(ac))) == -2) + 1; lp[which.max(ac[lp])]
  } else NA
  rs <- approx(ti, fi, xout = seq(a, b, length.out = n_resample))$y
  c(
    dur = dur, n_pts = sum(inside),
    f_mean = mean(fi), f_sd = sd(fi), f_min = min(fi), f_max = max(fi), f_range = diff(range(fi)),
    f_start = fi[1], f_end = fi[n], f_delta = fi[n] - fi[1],
    slope = unname(coef(lin)[2]), slope_rel = unname(coef(lin)[2]) * dur,
    curv = unname(quad), resid_sd = sd(res), resid_range = diff(range(res)),
    n_extrema = n_ext, osc_rate = n_ext / 2 / dur,
    period_ms = period * 1000, n_peaks = length(pk), amp_osc = if (length(pk) && length(tr)) mean(fd[pk]) - mean(fd[tr]) else 0,
    ac_period_ms = if (is.na(ac_lag)) NA else ac_lag * (ti[2] - ti[1]) * 1000,
    ac_peak = if (is.na(ac_lag)) NA else ac[ac_lag],
    tv_rate = sum(abs(d1)) / dur, tv_rel = sum(abs(d1)) / max(diff(range(fi)), 1e-3),
    setNames(rs - mean(rs), paste0("shape", seq_len(n_resample)))
  )
}

rows <- lapply(split(segs, segs$ID), function(s) {
  ct <- contours[contours$song.type == s$ID[1], ]
  if (nrow(ct) != 1) return(NULL)
  ff <- unlist(ct[, ff_cols])
  tt <- seq(ct$start, ct$end, length.out = length(ff))
  feats <- t(mapply(function(a, b) contour_features(tt, ff, a, b), s$start, s$end))
  data.frame(
    ID = s$ID, lek.song.type.year = s$lek.song.type.year, segment.pos = seq_len(nrow(s)),
    n_segments = nrow(s), rel_pos = (s$start - ct$start) / (ct$end - ct$start),
    song_f_mean = mean(ff), Category = s$Category, start = s$start, end = s$end, feats
  )
})
X <- do.call(rbind, rows)

# lek-song type grouping for cross-validation (e.g. BR1-A1 from ID 1-BR1-A1-2011)
X$lek.song.type <- sub("-[0-9]{4}$", "", sub("^[0-9]+-", "", X$ID))
X$f_mean_rel <- X$f_mean - X$song_f_mean

# neighbour context: previous and next segment within the same song
ctx <- c("slope", "slope_rel", "f_mean", "f_range", "osc_rate", "resid_sd", "dur", "tv_rate")
for (k in c(-1, 1)) {
  nb <- do.call(rbind, lapply(split(X, X$ID), function(s) {
    idx <- seq_len(nrow(s)) + k
    m <- s[ifelse(idx < 1 | idx > nrow(s), NA, idx), ctx]
    m$f_mean <- m$f_mean - s$f_mean
    m
  }))
  names(nb) <- paste0(ifelse(k < 0, "prev_", "next_"), ctx)
  X <- cbind(X, nb[match(rownames(X), rownames(X)), , drop = FALSE])
}

dir.create("./data/processed/segment_classifier", showWarnings = FALSE)
write.csv(X, "./data/processed/segment_classifier/segment_features.csv", row.names = FALSE)
cat(nrow(X), "segments from", length(unique(X$ID)), "songs and", length(unique(X$lek.song.type)), "song types\n")
print(table(X$Category))
