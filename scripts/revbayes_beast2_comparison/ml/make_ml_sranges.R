# Write the BEAST2 marginal likelihood runs of the stratigraphic range model (sRanges) for the 5 ms pilot,
# to compare with the sampled ancestor (SA) runs of make_ml.R (bounds, relaxed clock). HC1 and BR1, 3 replicates.
#
# Two versions:
#   SR      : first and last year of a song type are the ends of one range (as in the pilot sRanges runs)
#   SRsingle: every sample is its own single-occurrence range. This is the SA model written as an sRanges
#             model: the sRanges density is defined on oriented trees and equals the SA density minus
#             ln(2) per bifurcation, which integrates back to the SA marginal likelihood once the
#             orientations are summed over. So SRsingle should match SA; any gap shows a normalisation
#             or sampling problem that would also bias SR vs SA.
#
# Both versions add two tree moves that the pilot sRanges XMLs did not have: LeftRightChildSwap (changes the
# orientation, which is part of the sRanges tree) and SRLeafToSampledAncestorJump (moves a single-occurrence
# range between tip and sampled ancestor, the sRanges version of the SA LeafToSampledAncestorJump).
# Everything else (continuous relaxed clock rates with an Exponential(1) prior, 100 powers, 1.2 million
# proposals per power) is as in make_ml.R.
#
# Run from the project root:
#   Rscript scripts/revbayes_beast2_comparison/ml/make_ml_sranges.R

leks <- c("HC1", "BR1")
versions <- c("SR", "SRsingle")
n_rep <- 3
ml_dir <- "./scripts/revbayes_beast2_comparison/ml"
work_b2 <- file.path(Sys.getenv("HOME"), "lbh_ml_work")
step_length <- 1200000
pre_burnin <- 1500000

for (l in leks) {
  pilot_xml <- readLines(file.path("./scripts/revbayes_beast2_comparison/beast2", paste0(l, "_SR.xml")))
  taxa <- sub('.*<taxon id="([^"]+)".*', "\\1", grep("<taxon id=", pilot_xml, value = TRUE))
  n_branches <- 2 * length(taxa) - 2

  for (v in versions) {
    id <- paste(l, "bounds", v, "relaxed", sep = "_")
    x <- pilot_xml

    if (v == "SRsingle") {
      # every sample its own range, in the tree and in the starting tree
      def <- grep('spec="sr.evolution.sranges.StratigraphicRange"', x)
      x <- c(x[seq_len(min(def) - 1)],
             sprintf('          <stratigraphicRange id="range_%s" spec="sr.evolution.sranges.StratigraphicRange" firstOccurrence="@%s" lastOccurrence="@%s"/>', taxa, taxa, taxa),
             x[(max(def) + 1):length(x)])
      ref <- grep('<stratigraphicRange idref=', x)
      x <- c(x[seq_len(min(ref) - 1)],
             sprintf('      <stratigraphicRange idref="range_%s"/>', taxa),
             x[(max(ref) + 1):length(x)])
    }

    # relaxed clock with continuous rates and their prior
    cat_line <- grep('id="rateCategories"', x)
    x[cat_line] <- sprintf('      <parameter id="rates" name="stateNode" lower="0.0" dimension="%d">1.0</parameter>', n_branches)
    x <- sub('rateCategories="@rateCategories"', 'rates="@rates"', x, fixed = TRUE)
    cm <- grep('x="@clockMean"', x)
    x <- append(x, '        <distribution spec="beast.base.inference.distribution.Prior" x="@rates"><distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/></distribution>', after = cm)
    x <- x[!grepl("@rateCategories", x)]
    op <- max(grep("DeltaExchangeOperator", x))
    x <- append(x, c('    <operator spec="ScaleOperator" parameter="@rates" scaleFactor="0.75" weight="30"/>',
                     '    <operator spec="UpDownOperator" scaleFactor="0.75" weight="5" up="@clockMean" down="@rates"/>'), after = op)

    # extra tree moves
    wb <- grep("SRWilsonBalding", x)
    x <- append(x, c('    <operator spec="sr.evolution.operators.LeftRightChildSwap" tree="@tree" weight="10"/>',
                     '    <operator spec="sr.evolution.operators.SRLeafToSampledAncestorJump" tree="@tree" weight="10"/>'), after = wb)

    for (r in seq_len(n_rep)) {
      y <- x
      l1 <- grep('<logger id="tracelog"', y); l2 <- grep("</run>", y)
      y <- c(y[seq_len(l1 - 1)],
             sprintf('    <logger id="tracelog" fileName="trace.log" logEvery="%d" model="@posterior" sanitiseHeaders="true">', step_length / 100),
             '      <log idref="posterior"/>', '      <log idref="likelihood"/>', '      <log idref="prior"/>',
             '      <log idref="clockMean"/>', '      <log idref="origin"/>',
             '    </logger>',
             '  </mcmc>',
             '</run>',
             y[(l2 + 1):length(y)])
      run_line <- grep('<run id="mcmc"', y)
      root <- file.path(work_b2, paste0(id, "_r", r))
      y <- c(y[seq_len(run_line - 1)],
             sprintf('  <run spec="modelselection.inference.PathSampler" nrOfSteps="100" alpha="0.2" chainLength="%.0f" burnInPercentage="60" preBurnin="%.0f" rootdir="%s" deleteOldLogs="true">',
                     step_length, pre_burnin, root),
             'cd $(dir)',
             '$(java) -cp $(java.class.path) beast.pkgmgmt.launcher.BeastLauncher $(resume/overwrite) -java -seed $(seed) beast.xml',
             sub('<run id="mcmc" spec="MCMC" chainLength="[0-9.e+]*"', sprintf('<mcmc id="mcmc" spec="MCMC" chainLength="%.0f"', step_length), y[run_line]),
             y[(run_line + 1):length(y)])
      writeLines(y, file.path(ml_dir, "beast2", paste0(id, "_r", r, ".xml")))
    }
    cat("wrote", id, "\n")
  }
}
