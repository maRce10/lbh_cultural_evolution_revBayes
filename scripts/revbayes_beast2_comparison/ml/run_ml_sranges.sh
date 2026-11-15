#!/bin/bash
# Run the sRanges marginal likelihood runs of make_ml_sranges.R (12 runs, 1 core each, all at once).
# Run from the project root:
#   bash scripts/revbayes_beast2_comparison/ml/run_ml_sranges.sh
beast="$HOME/beast2_7/beast/bin/beast"
out="output/revbayes_beast2_comparison/ml"
mkdir -p "$out/beast2"
[ -f "$out/run_times.csv" ] || echo "program,model,exit_code,seconds" > "$out/run_times.csv"
run_b2() {
  id="$1"; seed=$(( 100 + ${id##*_r} )); start=$(date +%s)
  "$beast" -seed "$seed" -overwrite "scripts/revbayes_beast2_comparison/ml/beast2/${id}.xml" > "$out/beast2/${id}_screen.txt" 2>&1
  echo "beast2,${id},$?,$(( $(date +%s) - start ))" >> "$out/run_times.csv"
}
export -f run_b2
export beast out
ids=$(for l in BR1 HC1; do for v in SR SRsingle; do for r in 1 2 3; do echo "${l}_bounds_${v}_relaxed_r${r}"; done; done; done)
echo "$ids" | xargs -P 12 -I{} bash -c 'run_b2 {}'
echo "done"
