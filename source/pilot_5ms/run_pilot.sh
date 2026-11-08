#!/bin/bash
# Run the 6 pilot models, 4 at a time, each with 8 processes (2 runs x 4 chains).
# Run from the project root, after source/pilot_5ms/setup_pilot.R:
#   bash source/pilot_5ms/run_pilot.sh [leks]
# RevBayes 1.4.2 from bioconda: mamba create -n revbayes -c conda-forge -c bioconda revbayes=1.4.2

B="$HOME/mambaforge/envs/revbayes/bin"
mkdir -p output/pilot_5ms/screen

run_model() {
  id="$1"
  start=$(date +%s)
  "$B/mpirun" --mca btl_vader_single_copy_mechanism none -np 8 "$B/rb-mpi" "source/pilot_5ms/run_${id}.Rev" \
    > "output/pilot_5ms/screen/${id}.txt" 2>&1
  status=$?
  echo "${id},${status},$(( $(date +%s) - start ))" >> output/pilot_5ms/run_times.csv
}
export -f run_model
export B

[ -f output/pilot_5ms/run_times.csv ] || echo "model,exit_code,seconds" > output/pilot_5ms/run_times.csv

# leks given as arguments (default: all three), smallest datasets first
for lek in "${@:-HC1}" $( [ $# -eq 0 ] && echo BR1 SUR ); do
  for set in bounds all; do
    echo "${lek}_5ms_${set}"
  done
done | xargs -P 4 -I{} bash -c 'run_model {}'

echo "done"
cat output/pilot_5ms/run_times.csv
