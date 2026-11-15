# Fine-resolution trill features from the audio, guided by the tailored contours
# Adds to segment_features.csv:
#  - a peak-frequency track at 0.25 ms hop, searched only within +/-1.5 kHz of the
#    (interpolated) tailored contour, so it cannot jump to other harmonics
#  - trill period from that track (peak-to-peak intervals and autocorrelation)
#  - amplitude-envelope modulation period and depth
# Output: data/processed/segment_classifier/segment_features_audio.csv

library(tuneR)
library(seewave)

X <- read.csv("./data/processed/segment_classifier/segment_features.csv")
contours <- read.csv("./data/raw/segment_classification/seltailor_output.csv")
ff_cols <- grep("^ff", names(contours))
wav <- readWave("./data/raw/segment_classification/contour-degradation_325_songs.wav")
sr <- wav@samp.rate
y <- ffilter(wav, f = sr, from = 1000, to = 12000, bandpass = TRUE, output = "matrix")[, 1]

WL <- 128; HOP <- 11; NFFT <- 1024             # 2.9 ms window, 0.25 ms hop, zero-padded
win <- 0.5 - 0.5 * cos(2 * pi * (0:(WL - 1)) / (WL - 1))
fbin <- (0:(NFFT / 2)) * sr / NFFT / 1000      # kHz

# median interval between successive local maxima above a minimum prominence
peak_period <- function(x, t, prom) {
  d <- x - stats::lowess(t, x, f = 0.5)$y
  pk <- which(diff(sign(diff(d))) == -2) + 1; pk <- pk[d[pk] > prom]
  tr <- which(diff(sign(diff(d))) == 2) + 1;  tr <- tr[d[tr] < -prom]
  ivl <- c(diff(t[pk]), diff(t[tr]))
  c(period = if (length(ivl)) median(ivl) * 1000 else NA, n = length(pk))
}

# first autocorrelation peak between 3 and 40 ms
ac_period <- function(x, dt) {
  if (length(x) < 20 || sd(x) == 0) return(c(NA, NA))
  x <- x - stats::lowess(seq_along(x), x, f = 0.5)$y
  ac <- acf(x, lag.max = min(length(x) - 2, round(0.04 / dt)), plot = FALSE)$acf[-1]
  lags <- seq_along(ac) * dt * 1000
  lp <- which(diff(sign(diff(ac))) == -2) + 1; lp <- lp[lags[lp] >= 3]
  if (!length(lp)) return(c(NA, NA))
  b <- lp[which.max(ac[lp])]; c(lags[b], ac[b])
}

feats <- t(sapply(seq_len(nrow(X)), function(i) {
  ct <- contours[contours$song.type == X$ID[i], ]
  ff <- unlist(ct[, ff_cols]); tt <- seq(ct$start, ct$end, length.out = length(ff))
  a <- X$start[i]; b <- X$end[i]
  starts <- seq(round(a * sr) - WL / 2, round(b * sr) - WL / 2, by = HOP)
  starts <- starts[starts > 0 & starts + WL <= length(y)]
  if (length(starts) < 4) return(rep(NA, 8))
  tc <- (starts + WL / 2) / sr
  guide <- approx(tt, ff, xout = tc, rule = 2)$y
  pf <- env <- numeric(length(starts))
  for (j in seq_along(starts)) {
    fr <- y[starts[j] + 0:(WL - 1)] * win
    sp <- Mod(fft(c(fr, rep(0, NFFT - WL))))[1:(NFFT / 2 + 1)]
    ok <- abs(fbin - guide[j]) <= 1.5
    pf[j] <- fbin[ok][which.max(sp[ok])]
    env[j] <- sqrt(mean(fr^2))
  }
  lenv <- log(env + 1e-9)
  pp <- peak_period(pf, tc, prom = 0.1)
  pa <- ac_period(pf, HOP / sr)
  ea <- ac_period(lenv, HOP / sr)
  c(af_period_ms = pp[["period"]], af_n_peaks = pp[["n"]], af_ac_period_ms = pa[1], af_ac_peak = pa[2],
    af_range = diff(range(pf)), env_ac_period_ms = ea[1], env_ac_peak = ea[2], env_depth = sd(lenv))
}))

X <- cbind(X, feats)
write.csv(X, "./data/processed/segment_classifier/segment_features_audio.csv", row.names = FALSE)
print(aggregate(cbind(af_period_ms, af_ac_period_ms, env_ac_period_ms) ~ Category,
                X[X$Category %in% c("F", "M", "S"), ], median))
