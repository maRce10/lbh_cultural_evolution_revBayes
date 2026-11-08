# Write BEAST2 XML for one lek of the pilot, with the bounds data and the same model as the RevBayes pilot
# (relaxed clock, 5-state GTR + G, same priors, rho = 1, conditioned on the origin time) under two tree priors:
#   SA: fossilized birth-death with sampled ancestors, where first and last year are separate samples
#   SR: stratigraphic ranges (sRanges), where first and last year mark the ends of one range
# The sequence of a song type is kept at both ends of its range, as in the sRanges examples.
#
# Run from the project root:
#   Rscript source/pilot_5ms/beast2/make_xml.R <lek> <chain length> <log every> [suffix] [dataset: bounds or all]
# With all, only the SA model is written (ranges only use the first and last year)

args <- commandArgs(TRUE)
lek <- args[1]
chain <- as.numeric(args[2])
log_every <- as.numeric(args[3])
suffix <- if (length(args) > 3) args[4] else ""
dataset <- if (length(args) > 4) args[5] else "bounds"

out_dir <- "./source/pilot_5ms/beast2"
log_dir <- "output/pilot_5ms/beast2"

nex <- readLines(file.path("./data/processed/pilot_5ms/nexus", paste0(lek, "_5ms_", dataset, ".nex")))
rows <- strsplit(grep("\t", nex, value = TRUE), "\t")
seqs <- setNames(sapply(rows, `[`, 2), sapply(rows, `[`, 1))
tax <- read.delim(file.path("./data/processed/pilot_5ms/taxa", paste0(lek, "_5ms_", dataset, ".tsv")))
ages <- setNames(tax$age, tax$taxon)[names(seqs)]
type <- sub("-[0-9]{4}$", "", names(seqs))
year <- as.numeric(sub("^.*-", "", names(seqs)))

taxa_xml <- paste0('      <taxon id="', names(seqs), '" spec="Taxon"/>', collapse = "\n")
dates <- paste0(names(seqs), "=", ages, collapse = ",")
seq_xml <- paste0('    <sequence taxon="', names(seqs), '" totalcount="5" value="', seqs, '"/>', collapse = "\n")
n_branches <- 2 * length(seqs) - 2

ranges <- do.call(rbind, lapply(split(names(seqs), type), function(x) {
  y <- as.numeric(sub("^.*-", "", x))
  data.frame(type = sub("-[0-9]{4}$", "", x[1]), first = x[which.min(y)], last = x[which.max(y)])
}))
range_xml <- paste0('          <stratigraphicRange id="range_', ranges$type, '" spec="sr.evolution.sranges.StratigraphicRange" firstOccurrence="@',
                    ranges$first, '" lastOccurrence="@', ranges$last, '"/>', collapse = "\n")
range_ref <- paste0('      <stratigraphicRange idref="range_', ranges$type, '"/>', collapse = "\n")

write_xml <- function(prior) {
  id <- paste0(lek, "_", prior, suffix)
  sr <- prior == "SR"

  tree <- if (sr) sprintf('      <tree id="tree" spec="sr.evolution.tree.SRTree" nodetype="sr.evolution.tree.SRNode" name="stateNode">
        <trait id="dates" spec="beast.base.evolution.tree.TraitSet" traitname="date-backward" value="%s">
          <taxa idref="taxonSet"/>
        </trait>
        <taxonset idref="taxonSet"/>
%s
      </tree>', dates, range_xml) else sprintf('      <tree id="tree" spec="beast.base.evolution.tree.Tree" name="stateNode">
        <trait id="dates" spec="beast.base.evolution.tree.TraitSet" traitname="date-backward" value="%s">
          <taxa idref="taxonSet"/>
        </trait>
        <taxonset idref="taxonSet"/>
      </tree>', dates)

  init <- if (sr) sprintf('    <init estimate="false" id="startTree" initial="@tree" spec="sr.evolution.tree.RandomSRangeTree" nodetype="sr.evolution.tree.SRNode" taxonset="@taxonSet">
      <populationModel spec="beast.base.evolution.tree.coalescent.ConstantPopulation"><parameter name="popSize">1.0</parameter></populationModel>
%s
    </init>', range_ref) else '    <init estimate="false" id="startTree" initial="@tree" spec="beast.base.evolution.tree.coalescent.RandomTree" taxonset="@taxonSet">
      <populationModel spec="beast.base.evolution.tree.coalescent.ConstantPopulation"><parameter name="popSize">1.0</parameter></populationModel>
    </init>'

  prior_spec <- if (sr) "sr.speciation.SRangesBirthDeathModel" else "sa.evolution.speciation.SABirthDeathModel"

  tree_ops <- if (sr) '    <operator spec="sr.evolution.operators.SRWilsonBalding" tree="@tree" weight="20"/>
    <operator spec="sa.evolution.operators.SAUniform" tree="@tree" weight="20"/>
    <operator spec="sa.evolution.operators.SAScaleOperator" scaleFactor="0.95" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAScaleOperator" scaleFactor="0.95" rootOnly="true" tree="@tree" weight="5"/>' else '    <operator spec="sa.evolution.operators.LeafToSampledAncestorJump" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAWilsonBalding" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAExchange" isNarrow="false" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAExchange" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAUniform" tree="@tree" weight="20"/>
    <operator spec="sa.evolution.operators.SAScaleOperator" scaleFactor="0.95" tree="@tree" weight="10"/>
    <operator spec="sa.evolution.operators.SAScaleOperator" scaleFactor="0.95" rootOnly="true" tree="@tree" weight="5"/>'

  sa_log <- if (sr) "" else '      <log spec="sa.evolution.tree.SampledAncestorLogger" tree="@tree"/>'

  xml <- sprintf('<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<beast version="2.7" namespace="beast.base.core:beast.base.inference:beast.base.evolution.alignment:beast.base.evolution.tree:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood:beast.base.inference.operator:beast.base.evolution.operator">

  <data id="song" spec="beast.base.evolution.alignment.Alignment" dataType="user defined">
    <userDataType spec="beast.base.evolution.datatype.UserDataType" codelength="1" states="5" codeMap="D=0,F=1,P=2,S=3,U=4,-=0 1 2 3 4,?=0 1 2 3 4"/>
%s
  </data>

  <taxonset id="taxonSet" spec="beast.base.evolution.alignment.TaxonSet">
%s
  </taxonset>

  <run id="mcmc" spec="MCMC" chainLength="%.0f">
    <state id="state" storeEvery="%.0f">
%s
      <parameter id="origin" name="stateNode" lower="0.0" upper="1000.0">50.0</parameter>
      <parameter id="birth" name="stateNode" lower="0.0">0.1</parameter>
      <parameter id="death" name="stateNode" lower="0.0">0.05</parameter>
      <parameter id="psi" name="stateNode" lower="0.0">0.1</parameter>
      <parameter id="clockMean" name="stateNode" lower="0.0">0.1</parameter>
      <stateNode id="rateCategories" spec="parameter.IntegerParameter" dimension="%d">1</stateNode>
      <parameter id="alpha" name="stateNode" lower="0.0">1.0</parameter>
      <parameter id="er" name="stateNode" lower="0.0" upper="1.0" dimension="10">0.1</parameter>
      <parameter id="sf" name="stateNode" lower="0.0" upper="1.0" dimension="5">0.2</parameter>
    </state>

%s

    <distribution id="posterior" spec="CompoundDistribution">
      <distribution id="prior" spec="CompoundDistribution">
        <distribution id="treePrior" spec="%s" tree="@tree" origin="@origin" birthRate="@birth" deathRate="@death" samplingRate="@psi" conditionOnSampling="false" conditionOnRhoSampling="false">
          <parameter estimate="false" name="removalProbability" lower="0.0" upper="1.0">0.0</parameter>
          <parameter estimate="false" name="rho" lower="0.0" upper="1.0">1.0</parameter>
        </distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@origin"><distr spec="beast.base.inference.distribution.Uniform" lower="0.0" upper="1000.0"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@birth"><distr spec="beast.base.inference.distribution.Exponential" mean="0.1"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@death"><distr spec="beast.base.inference.distribution.Exponential" mean="0.1"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@psi"><distr spec="beast.base.inference.distribution.Exponential" mean="0.1"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@clockMean"><distr spec="beast.base.inference.distribution.Exponential" mean="0.1"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@alpha"><distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@er"><distr spec="beast.base.inference.distribution.Dirichlet" alpha="1 1 1 1 1 1 1 1 1 1"/></distribution>
        <distribution spec="beast.base.inference.distribution.Prior" x="@sf"><distr spec="beast.base.inference.distribution.Dirichlet" alpha="1 1 1 1 1"/></distribution>
      </distribution>
      <distribution id="likelihood" spec="TreeLikelihood" data="@song" tree="@tree">
        <siteModel id="siteModel" spec="SiteModel" gammaCategoryCount="4" shape="@alpha">
          <substModel id="gtr" spec="beastclassic.evolution.substitutionmodel.SVSGeneralSubstitutionModel" symmetric="true" rates="@er">
            <rateIndicator spec="parameter.BooleanParameter" estimate="false" dimension="10">true</rateIndicator>
            <frequencies spec="Frequencies" frequencies="@sf"/>
          </substModel>
        </siteModel>
        <branchRateModel id="clock" spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel" clock.rate="@clockMean" rateCategories="@rateCategories" tree="@tree">
          <distr spec="beast.base.inference.distribution.Exponential" mean="1.0"/>
        </branchRateModel>
      </distribution>
    </distribution>

%s
    <operator spec="ScaleOperator" parameter="@origin" scaleFactor="0.75" weight="5"/>
    <operator spec="ScaleOperator" parameter="@birth" scaleFactor="0.75" weight="5"/>
    <operator spec="ScaleOperator" parameter="@death" scaleFactor="0.75" weight="5"/>
    <operator spec="ScaleOperator" parameter="@psi" scaleFactor="0.75" weight="5"/>
    <operator spec="ScaleOperator" parameter="@clockMean" scaleFactor="0.75" weight="5"/>
    <operator spec="ScaleOperator" parameter="@alpha" scaleFactor="0.75" weight="5"/>
    <operator spec="DeltaExchangeOperator" parameter="@er" delta="0.05" weight="20"/>
    <operator spec="DeltaExchangeOperator" parameter="@sf" delta="0.05" weight="20"/>
    <operator spec="IntRandomWalkOperator" parameter="@rateCategories" windowSize="1" weight="10"/>
    <operator spec="SwapOperator" intparameter="@rateCategories" weight="10"/>
    <operator spec="UniformOperator" parameter="@rateCategories" weight="10"/>

    <logger id="tracelog" fileName="%s/%s.$(seed).log" logEvery="%.0f" model="@posterior" sanitiseHeaders="true" sort="smart">
      <log idref="posterior"/>
      <log idref="likelihood"/>
      <log idref="prior"/>
      <log idref="origin"/>
      <log idref="birth"/>
      <log idref="death"/>
      <log idref="psi"/>
      <log idref="clockMean"/>
      <log idref="alpha"/>
      <log idref="er"/>
      <log idref="sf"/>
      <log spec="beast.base.evolution.tree.TreeStatLogger" tree="@tree"/>
%s
    </logger>
    <logger id="treelog" fileName="%s/%s.$(seed).trees" logEvery="%.0f" mode="tree">
      <log spec="beast.base.evolution.TreeWithMetaDataLogger" tree="@tree" branchratemodel="@clock"/>
    </logger>
    <logger id="screenlog" logEvery="%.0f"><log idref="posterior"/><log idref="likelihood"/></logger>
  </run>
</beast>
', seq_xml, taxa_xml, chain, log_every, tree, n_branches, init, prior_spec, tree_ops,
     log_dir, id, log_every, sa_log, log_dir, id, log_every, max(log_every, chain / 20))
  writeLines(xml, file.path(out_dir, paste0(id, ".xml")))
  cat("wrote", id, "\n")
}

write_xml("SA")
if (dataset == "bounds") write_xml("SR")
