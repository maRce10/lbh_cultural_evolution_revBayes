#!/bin/bash
# =============================================================================
# PPS.sh
#
# Posterior predictive simulation (absolute model fit) under the relaxed clock
# (Uexp) for every lek x alignment x historical sampling scheme x historical
# record combination with an alignment and a taxa file.
#
# For each dataset, three stages, each a separate RevBayes run:
#   1. mcmc   posterior MCMC                 (pps_model_Uexp.Rev + pps_MCMC.Rev)
#   2. sim    simulate data from posterior   (pps_model_Uexp.Rev + pps_Simulation.Rev)
#   3. stats  summary statistics, p-values   (pps_DataSummary.Rev + pps_PValues.Rev)
# No dataset-specific Rev files are needed: the dataset's file names and extant
# clade are passed to RevBayes as variables, and the full command is written
# at the top of each log.
#
# OUTPUT, one folder per dataset: $OUT_DIR/<dataset>/
#   <dataset>.log/.trees/_posterior.var    MCMC output
#   sims/posterior_predictive_sim_<i>/     simulated data
#   simulated_data.csv, empirical_data.csv summary statistics
#   data_pvalues_effectsizes.csv           p-values and effect sizes
#   logs/<stage>.log                       RevBayes output of each stage
# and across datasets:
#   $OUT_DIR/pps_results.tsv   p-values and effect sizes, one row per dataset
#                              and statistic (rebuilt at the end of every run)
#   $OUT_DIR/pps_timing.tsv    wall time of every stage of every dataset
#   $OUT_DIR/missing.txt       combinations skipped, with the reason
#
# USAGE (run from the project home directory, the one containing scripts/ and data/):
#   bash scripts/PPS.sh                       # all stages, skipping finished ones
#   STAGES="mcmc" bash scripts/PPS.sh         # only some stages
#   STAGES="stats" bash scripts/PPS.sh        # e.g. recompute statistics only
#   ONLY="BR1_agnostic" bash scripts/PPS.sh   # only datasets whose name contains this
#   FORCE=1 bash scripts/PPS.sh               # rerun stages even if finished
# =============================================================================

set -uo pipefail

# ----------------------------- configuration --------------------------------
RB_DIR="/programs/revbayes-v1.4.2"
#RB_DIR="/programs/revbayes-mpi-1.2.1"
RB_EXE="${RB_DIR}/rb-mpi"
NPROC=8                                     # MPI processes per RevBayes run
SCRIPT_DIR="./scripts"
ALIGN_DIR="./data/processed/nexus"
TAXA_DIR="./data/processed/taxa"
OUT_DIR="./output/revbayes/pps"
STAGES="${STAGES:-mcmc sim stats}"
ONLY="${ONLY:-}"
FORCE="${FORCE:-0}"

leks=("BR1" "CCE" "HC1" "SUR" "TR1")
aligns=("agnostic" "optimal" "prank")
historical=("recent" "full" "none")
times=("all" "bounds" "extant")
CLOCK="Uexp"
SEED=1101

# MCMC length per lek (generations). The old MCMC.Rev used these values for
# SUR and CCE; shorter runs may be enough for the smaller leks.
declare -A BURNIN_GENS=( ["BR1"]=100000 ["CCE"]=100000 ["HC1"]=100000 ["SUR"]=100000 ["TR1"]=100000 )
declare -A RUN_GENS=(    ["BR1"]=300000 ["CCE"]=300000 ["HC1"]=300000 ["SUR"]=300000 ["TR1"]=300000 )
PRINTGEN=100        # sampling frequency: RUN_GENS/PRINTGEN posterior samples
PPS_THINNING=3      # simulate from every 3rd sample: 300000/100/3 = 1000 data sets

# Input file names (sampling times before historical record)
align_name () { echo "${lek}_${aln}_${tim}_${hist}_alignment.nex"; }
taxa_name ()  { echo "${lek}_${tim}_${hist}_taxa.tsv"; }

# ------------------------------- helpers ------------------------------------
fmt_time () { printf "%02d:%02d:%02d" $(( $1/3600 )) $(( $1%3600/60 )) $(( $1%60 )); }

# Extant taxa = taxa sampled in the last year (last 4 digits of the name)
extant_taxa () {
  awk -F'\t' '
    NR == 1 { next }
    { gsub(/\r/, ""); name = $1
      if (match(name, /[0-9][0-9][0-9][0-9]$/)) {
        yr = substr(name, RSTART, 4) + 0; taxon[NR] = name; year[NR] = yr
        if (yr > maxyr) maxyr = yr } }
    END { sep = ""
          for (i = 2; i <= NR; i++)
            if ((i in year) && year[i] == maxyr) { printf "%s\"%s\"", sep, taxon[i]; sep = ", " } }
  ' "$1"
}

# Run one stage: run_stage <stage> <rev commands>
run_stage () {
  local stage=$1 cmd=$2 log="${ddir}/logs/${1}.log" t0=$SECONDS
  echo "[$(date +%H:%M:%S)] ${dataset}: ${stage}"
  { echo "# COMMAND: ${cmd}"; echo; } > "$log"
  # q() makes RevBayes exit when the stage is done; under MPI it does not
  # reliably quit on its own when the piped commands run out
  { echo "$cmd"; echo "q()"; } | mpirun -np "$NPROC" "$RB_EXE" >> "$log" 2>&1
  local rb_exit=${PIPESTATUS[1]} secs=$(( SECONDS - t0 ))
  printf "%s\t%s\t%s\t%s\n" "$dataset" "$stage" "$secs" "$rb_exit" >> "$TIMING"
  echo "    finished in $(fmt_time $secs) (exit ${rb_exit}); log: ${log}"
  if grep -q "Error:" "$log"; then
    echo "    WARNING: RevBayes reported an error in ${log}" >&2; return 1
  fi
  return 0
}

wants () { [[ " $STAGES " == *" $1 "* ]]; }

# ------------------------------- setup --------------------------------------
[[ -x "$RB_EXE" ]] || { echo "ERROR: RevBayes not found at $RB_EXE" >&2; exit 1; }
mkdir -p "$OUT_DIR"
MISSING="${OUT_DIR}/missing.txt"; : > "$MISSING"
TIMING="${OUT_DIR}/pps_timing.tsv"
[[ -s "$TIMING" ]] || printf "dataset\tstage\tseconds\texit\n" > "$TIMING"

total_start=$SECONDS
n_run=0; n_skip=0; n_fail=0
echo "Started: $(date)   stages: ${STAGES}"

# ------------------------------- main loop ----------------------------------
for lek in "${leks[@]}"; do
 for aln in "${aligns[@]}"; do
  for hist in "${historical[@]}"; do
   for tim in "${times[@]}"; do

    dataset="${lek}_${aln}_${hist}_${tim}_${CLOCK}"
    [[ -z "$ONLY" || "$dataset" == *"$ONLY"* ]] || continue

    align_file="${ALIGN_DIR}/$(align_name)"
    taxa_file="${TAXA_DIR}/$(taxa_name)"
    if [[ ! -f "$align_file" ]]; then echo -e "${dataset}\tno alignment: ${align_file}" >> "$MISSING"; continue; fi
    if [[ ! -f "$taxa_file"  ]]; then echo -e "${dataset}\tno taxa file: ${taxa_file}"  >> "$MISSING"; continue; fi
    extant=$(extant_taxa "$taxa_file")
    if [[ -z "$extant" ]]; then echo -e "${dataset}\tno extant taxa in ${taxa_file}" >> "$MISSING"; continue; fi

    ddir="${OUT_DIR}/${dataset}"
    base="${ddir}/${dataset}"
    sim_dir="${ddir}/sims"
    mkdir -p "${ddir}/logs"

    # Variables shared by all stages
    vars="seed(${SEED});
alignment = \"${align_file}\"; fossils = \"${taxa_file}\";
base_str = \"${base}\"; output_file_name = base_str; state_name = base_str + \".state\";
sim_dir = \"${sim_dir}\";
clade_extant = clade(${extant});"

    # ---- stage 1: MCMC
    if wants mcmc; then
      if [[ -f "${base}.mcmc_done" && "$FORCE" != "1" ]]; then
        echo "[skip] ${dataset}: mcmc already finished"; n_skip=$(( n_skip + 1 ))
      else
        rm -f "${base}.mcmc_done" "${base}.sim_done"
        cmd="${vars}
printgen = ${PRINTGEN}; burnin_gens = ${BURNIN_GENS[$lek]}; run_gens = ${RUN_GENS[$lek]};
source(\"${SCRIPT_DIR}/pps_model_Uexp.Rev\");
source(\"${SCRIPT_DIR}/pps_MCMC.Rev\");"
        if run_stage mcmc "$cmd" && [[ -f "${base}.mcmc_done" ]]; then n_run=$(( n_run + 1 )); else n_fail=$(( n_fail + 1 )); continue; fi
      fi
    fi

    # ---- stage 2: simulation (needs a finished MCMC)
    if wants sim; then
      if [[ ! -f "${base}.mcmc_done" ]]; then
        echo "[skip] ${dataset}: sim needs a finished mcmc"; n_skip=$(( n_skip + 1 )); continue
      elif [[ -f "${base}.sim_done" && "$FORCE" != "1" ]]; then
        echo "[skip] ${dataset}: sim already finished"; n_skip=$(( n_skip + 1 ))
      else
        rm -rf "$sim_dir" "${base}.sim_done"
        cmd="${vars}
pps_thinning = ${PPS_THINNING};
source(\"${SCRIPT_DIR}/pps_model_Uexp.Rev\");
source(\"${SCRIPT_DIR}/pps_Simulation.Rev\");"
        if run_stage sim "$cmd" && [[ -f "${base}.sim_done" ]]; then n_run=$(( n_run + 1 )); else n_fail=$(( n_fail + 1 )); continue; fi
      fi
    fi

    # ---- stage 3: summary statistics and p-values (needs finished simulations)
    if wants stats; then
      if [[ ! -f "${base}.sim_done" ]]; then
        echo "[skip] ${dataset}: stats needs finished simulations"; n_skip=$(( n_skip + 1 )); continue
      fi
      n_sims=$(find "$sim_dir" -mindepth 1 -maxdepth 1 -type d -name 'posterior_predictive_sim_*' | wc -l)
      cmd="${vars}
num_post_sims = ${n_sims};
sim_pps_file = \"${ddir}/simulated_data.csv\";
emp_pps_file = \"${ddir}/empirical_data.csv\";
outfileName  = \"${ddir}/data_pvalues_effectsizes.csv\";
source(\"${SCRIPT_DIR}/pps_DataSummary.Rev\");
source(\"${SCRIPT_DIR}/pps_PValues.Rev\");"
      if run_stage stats "$cmd" && [[ -s "${ddir}/data_pvalues_effectsizes.csv" ]]; then
        n_run=$(( n_run + 1 ))
        echo "    ${n_sims} simulated data sets summarized"
      else
        n_fail=$(( n_fail + 1 ))
      fi
    fi

   done
  done
 done
done

# --------------------- combine p-values across datasets ---------------------
RESULTS="${OUT_DIR}/pps_results.tsv"
printf "lek\talignment\thistorical\ttimes\tclock\tstatistic\tlower_1tailed\tupper_1tailed\ttwo_tailed\tmidpoint\teffect_size\n" > "$RESULTS"
for f in "${OUT_DIR}"/*/data_pvalues_effectsizes.csv; do
  [[ -f "$f" ]] || continue
  ds=$(basename "$(dirname "$f")")
  IFS=_ read -r r_lek r_aln r_hist r_tim r_clk <<< "$ds"
  awk -F',' -v OFS='\t' -v pre="${r_lek}\t${r_aln}\t${r_hist}\t${r_tim}\t${r_clk}" '
    NR > 1 && NF >= 6 { gsub(/\r/, ""); print pre, $1, $2, $3, $4, $5, $6 }' "$f" >> "$RESULTS"
done

# ------------------------------- summary ------------------------------------
total=$(( SECONDS - total_start ))
echo "Finished: $(date)"
echo "Total wall time: $(fmt_time $total)"
echo "Stages run: ${n_run} | skipped: ${n_skip} | failed: ${n_fail}"
echo "Combinations skipped for missing inputs: $(wc -l < "$MISSING") (see $MISSING)"
awk -F'\t' 'NR>1 && $4==0 { n[$2]++; s[$2]+=$3 }
  END { for (k in n) printf "Mean time per dataset, %-5s: %s s over %d runs\n", k, int(s[k]/n[k]), n[k] }' "$TIMING"
echo "P-values and effect sizes for all datasets: $RESULTS"
