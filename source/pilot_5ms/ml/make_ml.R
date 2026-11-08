# Write the marginal likelihood runs of the 5 ms pilot: HC1 and BR1, all and bounds, global and relaxed clock,
# in RevBayes (power posterior, as in the first version) and BEAST2 (MODEL_SELECTION path sampler, SA tree prior).
# Both use 100 powers spaced as quantiles of Beta(0.2, 1) and report the stepping-stone estimate.
#
# BEAST2 relaxed clock: continuous branch rates with an Exponential(1) prior (times clockMean), as in RevBayes.
# The rate categories of the pilot XML have no prior term, so they would inflate the marginal likelihood
# by about n log(n) for n branches.
#
# BEAST2 effort per power is matched to RevBayes in number of proposals: RevBayes does about 300 proposals per
# generation and 4000 generations per stone (2000 pre-burnin + 2000), so 1.2 million proposals per BEAST2 step.
#
# Run from the project root:
#   Rscript source/pilot_5ms/ml/make_ml.R

leks <- c("HC1", "BR1")
datasets <- c("bounds", "all")
clocks <- c("global", "relaxed")
seeds <- c(101, 102, 103)

ml_dir <- "./source/pilot_5ms/ml"
out_rb <- "./output/pilot_5ms/ml/revbayes"
work_b2 <- file.path(Sys.getenv("HOME"), "lbh_ml_work")   # path sampler step folders, outside Dropbox
for (d in c(file.path(ml_dir, c("revbayes", "beast2")), out_rb, "./output/pilot_5ms/ml/beast2"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

clock_file <- c(global = "./source/model_Global.Rev", relaxed = "./source/model_UExp.Rev")
step_length <- 1200000
pre_burnin <- 1500000

for (l in leks) for (s in datasets) {
  pilot_rev <- readLines(file.path("./source/pilot_5ms", paste0("run_", l, "_5ms_", s, ".Rev")))
  head_lines <- grep("^(alignment|fossils|clade_extant) ", pilot_rev, value = TRUE)
  pilot_xml <- readLines(file.path("./source/pilot_5ms/beast2", paste0(l, "_SA", if (s == "all") "_all" else "", ".xml")))
  n_taxa <- length(grep("<taxon id=", pilot_xml))
  n_branches <- 2 * n_taxa - 2

  for (k in clocks) {
    id <- paste(l, s, k, sep = "_")

    # RevBayes, one file per replicate
    for (r in seq_along(seeds)) {
      writeLines(c(head_lines,
                   sprintf('base_str = "%s/%s_r%d"', out_rb, id, r),
                   sprintf('clock_model = "%s"', clock_file[[k]]),
                   sprintf("rep_seed = %d", seeds[r]),
                   'source("./source/pilot_5ms/ml/control_ml.Rev")'),
                 file.path(ml_dir, "revbayes", paste0(id, "_r", r, ".Rev")))
    }

    # BEAST2, one XML per replicate (each needs its own step folder)
    for (r in seq_along(seeds)) {
      x <- pilot_xml
      cat_line <- grep('id="rateCategories"', x)
      if (k == "relaxed") {
        x[cat_line] <- sprintf('      <parameter id="rates" name="stateNode" lower="0.0" dimension="%d">1.0</parameter>', n_branches)
        x <- sub('rateCategories="@rateCategories"', 'rates="@rates"', x, fixed = TRUE)
        cm <- grep('x="@clockMean"', x)
        x <- append(x, '        <distribution spec="beast.base.inference.distribution.Prior" x="@rates"><distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/></distribution>', after = cm)
      } else {
        x <- x[-cat_line]
        b1 <- grep("<branchRateModel", x); b2 <- grep("</branchRateModel>", x)
        x <- c(x[seq_len(b1 - 1)],
               '        <branchRateModel id="clock" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockMean"/>',
               x[(b2 + 1):length(x)])
      }
      # operators on rate categories
      x <- x[!grepl("@rateCategories", x)]
      if (k == "relaxed") {
        op <- max(grep("DeltaExchangeOperator", x))
        x <- append(x, c('    <operator spec="ScaleOperator" parameter="@rates" scaleFactor="0.75" weight="30"/>',
                         '    <operator spec="UpDownOperator" scaleFactor="0.75" weight="5" up="@clockMean" down="@rates"/>'), after = op)
      }
      # loggers: only a thinned trace log in each step; the path sampler adds its own likelihood log
      l1 <- grep('<logger id="tracelog"', x); l2 <- grep("</run>", x)
      x <- c(x[seq_len(l1 - 1)],
             sprintf('    <logger id="tracelog" fileName="trace.log" logEvery="%d" model="@posterior" sanitiseHeaders="true">', step_length / 100),
             '      <log idref="posterior"/>', '      <log idref="likelihood"/>', '      <log idref="prior"/>',
             '      <log idref="clockMean"/>', '      <log idref="origin"/>',
             '    </logger>',
             '  </mcmc>',
             '</run>',
             x[(l2 + 1):length(x)])
      run_line <- grep('<run id="mcmc"', x)
      root <- file.path(work_b2, paste0(id, "_r", r))
      x <- c(x[seq_len(run_line - 1)],
             sprintf('  <run spec="modelselection.inference.PathSampler" nrOfSteps="100" alpha="0.2" chainLength="%.0f" burnInPercentage="60" preBurnin="%.0f" rootdir="%s" deleteOldLogs="true">',
                     step_length, pre_burnin, root),
             'cd $(dir)',
             '$(java) -cp $(java.class.path) beast.pkgmgmt.launcher.BeastLauncher $(resume/overwrite) -java -seed $(seed) beast.xml',
             sub('<run id="mcmc" spec="MCMC" chainLength="[0-9.e+]*"', sprintf('<mcmc id="mcmc" spec="MCMC" chainLength="%.0f"', step_length), x[run_line]),
             x[(run_line + 1):length(x)])
      writeLines(x, file.path(ml_dir, "beast2", paste0(id, "_r", r, ".xml")))
    }
    cat("wrote", id, "\n")
  }
}
