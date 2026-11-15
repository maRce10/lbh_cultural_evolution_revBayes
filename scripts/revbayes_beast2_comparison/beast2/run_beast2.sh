#!/bin/bash
# Run one lek under the SA and sRanges tree priors, with two replicates each.
# Run from the project root, after scripts/revbayes_beast2_comparison/beast2/make_xml.R:
#   bash scripts/revbayes_beast2_comparison/beast2/run_beast2.sh <lek> [suffix, e.g. _all] [models, e.g. "SA"]
# BEAST 2.7.7 with the packages SA, MM, BEASTLabs, BEAST_CLASSIC and SRanges v0.1.3

lek="$1"
suffix="${2:-}"
models="${3:-SA SR}"
beast="$HOME/beast2_7/beast/bin/beast"
mkdir -p output/revbayes_beast2_comparison/beast2/screen
for m in $models; do
  for seed in 101 102; do
    "$beast" -seed "$seed" -overwrite -statefile "output/revbayes_beast2_comparison/beast2/screen/${lek}_${m}${suffix}_${seed}.state" "scripts/revbayes_beast2_comparison/beast2/${lek}_${m}${suffix}.xml" \
      > "output/revbayes_beast2_comparison/beast2/screen/${lek}_${m}${suffix}_${seed}.txt" 2>&1 &
  done
done
wait
echo finished
grep -l "Done!" output/revbayes_beast2_comparison/beast2/screen/${lek}_*${suffix}_*.txt
