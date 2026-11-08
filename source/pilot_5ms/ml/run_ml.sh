#!/bin/bash
# Run the marginal likelihood runs of the 5 ms pilot in RevBayes and BEAST2 at the same time:
# RevBayes 3 runs at a time with 8 processes each, BEAST2 8 runs at a time (1 core each).
# Run from the project root, after source/pilot_5ms/ml/make_ml.R:
#   bash source/pilot_5ms/ml/run_ml.sh

B="$HOME/mambaforge/envs/revbayes/bin"
beast="$HOME/beast2_7/beast/bin/beast"
out="output/pilot_5ms/ml"
mkdir -p "$out/revbayes" "$out/beast2"
[ -f "$out/run_times.csv" ] || echo "program,model,exit_code,seconds" > "$out/run_times.csv"

run_rb() {
  id="$1"; start=$(date +%s)
  "$B/mpirun" --mca btl_vader_single_copy_mechanism none -np 8 "$B/rb-mpi" "source/pilot_5ms/ml/revbayes/${id}.Rev" \
    > "$out/revbayes/${id}_screen.txt" 2>&1
  echo "revbayes,${id},$?,$(( $(date +%s) - start ))" >> "$out/run_times.csv"
}
run_b2() {
  id="$1"; seed=$(( 100 + ${id##*_r} )); start=$(date +%s)
  "$beast" -seed "$seed" -overwrite "source/pilot_5ms/ml/beast2/${id}.xml" > "$out/beast2/${id}_screen.txt" 2>&1
  echo "beast2,${id},$?,$(( $(date +%s) - start ))" >> "$out/run_times.csv"
}
export -f run_rb run_b2
export B beast out

ids=$(for l in BR1 HC1; do for s in all bounds; do for k in relaxed global; do for r in 1 2 3; do echo "${l}_${s}_${k}_r${r}"; done; done; done; done)

echo "$ids" | xargs -P 3 -I{} bash -c 'run_rb {}' &
echo "$ids" | xargs -P 8 -I{} bash -c 'run_b2 {}' &
wait
echo "done"
cat "$out/run_times.csv"
