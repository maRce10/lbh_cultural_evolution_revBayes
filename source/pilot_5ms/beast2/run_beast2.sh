#!/bin/bash
# Run one lek under the SA and sRanges tree priors, with two replicates each.
# Run from the project root, after source/pilot_5ms/beast2/make_xml.R:
#   bash source/pilot_5ms/beast2/run_beast2.sh <lek> [suffix, e.g. _all] [models, e.g. "SA"]
# BEAST 2.7.7 with the packages SA, MM, BEASTLabs, BEAST_CLASSIC and SRanges v0.1.3

lek="$1"
suffix="${2:-}"
models="${3:-SA SR}"
beast="$HOME/beast2_7/beast/bin/beast"
mkdir -p output/pilot_5ms/beast2/screen
for m in $models; do
  for seed in 101 102; do
    "$beast" -seed "$seed" -overwrite -statefile "output/pilot_5ms/beast2/screen/${lek}_${m}${suffix}_${seed}.state" "source/pilot_5ms/beast2/${lek}_${m}${suffix}.xml" \
      > "output/pilot_5ms/beast2/screen/${lek}_${m}${suffix}_${seed}.txt" 2>&1 &
  done
done
wait
echo finished
grep -l "Done!" output/pilot_5ms/beast2/screen/${lek}_*${suffix}_*.txt
