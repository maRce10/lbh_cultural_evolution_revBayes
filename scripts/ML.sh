#!/bin/bash
# =============================================================================
# ML.sh
#
# Approximates marginal likelihoods (path sampling and stepping-stone) in
# RevBayes for every combination of:
#     lek x alignment x historical record x sampling times x clock model x replicate
#
# STEP 1 - write the dataset-specific Rev scripts into ./scripts/lek_models_ML/,
#          one per combination whose alignment and taxa files both exist.
#          Each script records its input files and output name, the clade of
#          extant song types (taxa sampled in the last year), and then sources
#          the clock-specific ML control script.
# STEP 2 - run each script for three replicates, collecting the estimates.
#
# OUTPUT (all under $OUT_DIR):
#   ml_results.tsv  one row per run: dataset info, seed, both marginal
#                   likelihood estimates, wall time and status. Ready to
#                   read straight into R / Python as a table.
#   logs/           full RevBayes stdout for every run (<dataset>.<rep>.log),
#                   so estimates can always be re-extracted or checked.
#   missing.txt     combinations skipped, with the reason (missing input file,
#                   no extant taxa).
#   timing.txt      total wall time and a per-run timing summary.
#
# USAGE (run from the project home directory, the one containing scripts/ and data/):
#   bash ML.sh                     # write the Rev scripts, then run everything
#   GENERATE_ONLY=1 bash ML.sh     # only write the Rev scripts (to inspect them)
#   RESUME=1 bash ML.sh            # skip runs that already finished
#
# ASSUMPTIONS - check against your files:
#   1. Input file names (note the order of the parts differs between them):
#        alignment: ./data/processed/nexus/<lek>_<align>_<times>_<historical>_alignment.nex
#        taxa:      ./data/processed/taxa/<lek>_<times>_<historical>_taxa.tsv
#      (set by align_name() and taxa_name() in the configuration section)
#   2. Taxa files are tab-separated with a header line and the taxon name in
#      the first column; names end in the 4-digit sampling year (BR1-A1-2014).
#   3. The control scripts are ./scripts/Control_ML_global.Rev and
#      ./scripts/Control_ML_Uexp.Rev.
#   4. The ML control script prints its estimates on labelled lines:
#          print("SS_lnML", ss.marginal())
#          print("PS_lnML", ps.marginal())
# =============================================================================

set -uo pipefail   # error on unset variables; a pipeline fails if any part fails
                   # (no `set -e`: one failed run should not stop the whole loop)

# ----------------------------- configuration --------------------------------
RB_DIR="/programs/revbayes-1.3.0"
RB_EXE="${RB_DIR}/rb-mpi"
NPROC=32                                     # MPI processes per run
MODEL_DIR="./scripts/lek_models_ML"           # where the generated Rev scripts go
ALIGN_DIR="./data/processed/nexus"
TAXA_DIR="./data/processed/taxa"
CONTROL_DIR="./scripts"                       # holds Control_ML_<clock>.Rev
RB_OUT_DIR="./output/revbayes"               # RevBayes output (inside Rev scripts)
OUT_DIR="./output/revbayes/marginal_likelihoods"
RESUME="${RESUME:-0}"                        # 1 = skip runs with a finished log
GENERATE_ONLY="${GENERATE_ONLY:-0}"          # 1 = write Rev scripts, then stop

leks=("BR1" "CCE" "HC1" "SUR" "TR1")
aligns=("agnostic" "optimal" "prank")
historical=("recent" "full" "none")
times=("all" "bounds" "extant")
clocks=("global" "Uexp")
reps=("R1" "R2" "R3")

# Input file names. The alignments put the sampling times before the
# historical record, in both the alignment and the taxa file names.
# Change the order here if the naming changes again.
align_name () { echo "${lek}_${aln}_${tim}_${hist}_alignment.nex"; }   # e.g. BR1_agnostic_all_recent_alignment.nex
taxa_name ()  { echo "${lek}_${tim}_${hist}_taxa.tsv"; }               # e.g. BR1_all_recent_taxa.tsv

# One fixed seed per replicate, so every run can be reproduced exactly.
# The seed is set before the Rev script runs, so the random number `rv` in
# the script (used in output file names) is reproducible too.
declare -A REP_SEED=( ["R1"]=1101 ["R2"]=2202 ["R3"]=3303 )

# ------------------------------- setup --------------------------------------
LOG_DIR="${OUT_DIR}/logs"
RESULTS="${OUT_DIR}/ml_results.tsv"
MISSING="${OUT_DIR}/missing.txt"
TIMING="${OUT_DIR}/timing.txt"
mkdir -p "$LOG_DIR" "$MODEL_DIR"

# hh:mm:ss formatter for the timing report
fmt_time () { printf "%02d:%02d:%02d" $(( $1/3600 )) $(( $1%3600/60 )) $(( $1%60 )); }

# Pull a single number out of a line; empty if there is none
extract_number () { grep -Eo '[-+]?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?' <<< "$1" | tail -1; }

# Extant taxa = taxa sampled in the last year present in the taxa file, where
# the year is the last 4 digits of the taxon name. Prints them as a Rev
# argument list:  "BR1-A1-2014", "BR1-C1-2014", ...
extant_taxa () {
  awk -F'\t' '
    NR == 1 { next }                                    # skip header
    { gsub(/\r/, ""); name = $1
      if (match(name, /[0-9][0-9][0-9][0-9]$/)) {
        yr = substr(name, RSTART, 4) + 0
        taxon[NR] = name; year[NR] = yr
        if (yr > maxyr) maxyr = yr
      } }
    END { sep = ""
          for (i = 2; i <= NR; i++)
            if ((i in year) && year[i] == maxyr) { printf "%s\"%s\"", sep, taxon[i]; sep = ", " } }
  ' "$1"
}

# =============================================================================
# STEP 1: write the dataset-specific Rev scripts
# =============================================================================
: > "$MISSING"
n_written=0; n_missing=0

for clk in "${clocks[@]}"; do
  control="${CONTROL_DIR}/Control_ML_${clk}.Rev"
  [[ -f "$control" ]] || echo "WARNING: control script $control not found; ${clk} runs will fail" >&2
done

for lek in "${leks[@]}"; do
 for aln in "${aligns[@]}"; do
  for hist in "${historical[@]}"; do
   for tim in "${times[@]}"; do

     align_file="${ALIGN_DIR}/$(align_name)"
     taxa_file="${TAXA_DIR}/$(taxa_name)"

     # Not every combination exists: record why and move on
     if [[ ! -f "$align_file" ]]; then
       for clk in "${clocks[@]}"; do echo -e "${lek}_${aln}_${hist}_${tim}_${clk}\tno alignment: ${align_file}" >> "$MISSING"; done
       n_missing=$(( n_missing + ${#clocks[@]} )); continue
     fi
     if [[ ! -f "$taxa_file" ]]; then
       for clk in "${clocks[@]}"; do echo -e "${lek}_${aln}_${hist}_${tim}_${clk}\tno taxa file: ${taxa_file}" >> "$MISSING"; done
       n_missing=$(( n_missing + ${#clocks[@]} )); continue
     fi

     extant=$(extant_taxa "$taxa_file")
     if [[ -z "$extant" ]]; then
       for clk in "${clocks[@]}"; do echo -e "${lek}_${aln}_${hist}_${tim}_${clk}\tno extant taxa found in ${taxa_file}" >> "$MISSING"; done
       n_missing=$(( n_missing + ${#clocks[@]} )); continue
     fi

     for clk in "${clocks[@]}"; do
       dataset="${lek}_${aln}_${hist}_${tim}_${clk}"
       # Overwritten on every run, so the scripts always match the data files
       cat > "${MODEL_DIR}/${dataset}.Rev" <<EOF
# Generated by ML.sh on $(date +%Y-%m-%d) - edit ML.sh, not this file
alignment       = "${align_file}"
fossils = "${taxa_file}"
rv = floor(runif(1,100000,999999)[1])
output_file_name = "${RB_OUT_DIR}/${dataset}"
base_str = output_file_name + "_" + rv
state_name = output_file_name + ".state"


# Make clade of the extant songs #
# all song types sampled in the last year of ${taxa_file##*/}
clade_extant = clade(${extant})

source("${CONTROL_DIR}/Control_ML_${clk}.Rev")
EOF
       n_written=$(( n_written + 1 ))
     done

   done
  done
 done
done

echo "Wrote ${n_written} Rev scripts to ${MODEL_DIR}; ${n_missing} combinations skipped (see ${MISSING})"
[[ "$GENERATE_ONLY" == "1" ]] && exit 0

# =============================================================================
# STEP 2: run every generated script, three replicates each
# =============================================================================
[[ -x "$RB_EXE" ]] || { echo "ERROR: RevBayes not found at $RB_EXE" >&2; exit 1; }

# Write the table header only if the file is new (keeps rows when resuming)
if [[ ! -s "$RESULTS" ]]; then
  printf "lek\talignment\thistorical\ttimes\tclock\treplicate\tseed\tSS_lnML\tPS_lnML\tseconds\tstatus\n" > "$RESULTS"
fi

total_start=$SECONDS
n_done=0; n_failed=0; n_skipped=0

echo "Started: $(date)" | tee "$TIMING"

for lek in "${leks[@]}"; do
 for aln in "${aligns[@]}"; do
  for hist in "${historical[@]}"; do
   for tim in "${times[@]}"; do
    for clk in "${clocks[@]}"; do

      dataset="${lek}_${aln}_${hist}_${tim}_${clk}"
      rev_script="${MODEL_DIR}/${dataset}.Rev"
      [[ -f "$rev_script" ]] || continue      # skipped in step 1

      for rep in "${reps[@]}"; do
        seed=${REP_SEED[$rep]}
        log="${LOG_DIR}/${dataset}.${rep}.log"

        # Resume: skip runs whose log is marked as finished
        if [[ "$RESUME" == "1" && -s "$log" ]] && grep -q "^# STATUS: ok" "$log"; then
          echo "[skip] ${dataset} ${rep} (already finished)"
          n_skipped=$(( n_skipped + 1 ))
          continue
        fi

        echo "[$(date +%H:%M:%S)] starting ${dataset} ${rep} (seed ${seed})"
        run_start=$SECONDS

        # Commands piped to RevBayes: set the seed and replicate name, then run
        # the dataset script. stderr is kept in the log too, for debugging.
        rb_command="seed(${seed}); rep = \"${rep}\"; source(\"${rev_script}\");"
        echo "$rb_command" | mpirun -np "$NPROC" "$RB_EXE" > "$log" 2>&1
        rb_exit=${PIPESTATUS[1]}

        secs=$(( SECONDS - run_start ))

        # Read the estimates from the labelled lines printed by the control
        # script. Matching the labels (not "the last two lines") avoids picking
        # up numbers from file names such as "BR1" in RevBayes' messages.
        ss=$(extract_number "$(tr -d '\r' < "$log" | grep -E '^[[:space:]>]*SS_lnML' | tail -1 | sed 's/.*SS_lnML//')")
        ps=$(extract_number "$(tr -d '\r' < "$log" | grep -E '^[[:space:]>]*PS_lnML' | tail -1 | sed 's/.*PS_lnML//')")

        if [[ $rb_exit -eq 0 && -n "$ss" && -n "$ps" ]]; then
          status="ok"; n_done=$(( n_done + 1 ))
        else
          status="failed(exit=${rb_exit})"; n_failed=$(( n_failed + 1 ))
          ss=${ss:-NA}; ps=${ps:-NA}
          echo "  WARNING: ${dataset} ${rep} failed or gave no estimates; see $log" >&2
        fi
        echo "# STATUS: ${status}" >> "$log"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
          "$lek" "$aln" "$hist" "$tim" "$clk" "$rep" "$seed" \
          "$ss" "$ps" "$secs" "$status" >> "$RESULTS"

        echo "  done in $(fmt_time $secs): SS = ${ss}, PS = ${ps}"
      done

    done
   done
  done
 done
done

# ------------------------------- summary ------------------------------------
total=$(( SECONDS - total_start ))
{
  echo "Finished: $(date)"
  echo "Total wall time: $(fmt_time $total) (${total} s)"
  echo "Runs completed: ${n_done} | failed: ${n_failed} | skipped (resume): ${n_skipped}"
  echo "Combinations without input files: ${n_missing} (listed in $MISSING)"
  # Per-run timing from the results table (successful runs only)
  awk -F'\t' 'NR>1 && $11=="ok" { n++; s+=$10; if(min==""||$10<min)min=$10; if($10>max)max=$10 }
    END { if(n) printf "Per-run time (s): mean %.0f, min %d, max %d over %d runs\n", s/n, min, max, n }' "$RESULTS"
  echo "Per-run times for every dataset and replicate are in column 'seconds' of $RESULTS"
} | tee -a "$TIMING"
