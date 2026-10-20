#!/bin/bash
set -euo pipefail

SETTINGS="input/sys_settings.dat"
IS_SETT="input/ising_settings.dat"

OUTP_DIR="output"
RES_DIR="../nsl_exercises/sim_results/ex6"

rm -rf "${RES_DIR}"/*

run() {
  ./source/run
}

NPART=50
T_MIN=0.5
T_DELTA=0.1
TSTEPS=15
TEMPS=()
for i in $(seq 0 "${TSTEPS}"); do
  TEMPS+=($(awk -v T_MIN="${T_MIN}" -v T_DELTA="${T_DELTA}" -v i="${i}" 'BEGIN{printf "%.2f", T_MIN + i * T_DELTA'}))
done
H=0.02
ALGOS=("METROPOLIS" "GIBBS")
NSTEPS=10000
NBLOCKS=100
NSTEPS_EQUIL=2000
NBLOCKS_EQUIL=16

sed -i "s/^MODEL .*/MODEL ISING/" "${SETTINGS}"
sed -i "s/^RESTARTED .*/RESTARTED 0/" "${SETTINGS}"
sed -i "s/^NBLOCKS .*/NBLOCKS ${NBLOCKS}/" "${SETTINGS}"
sed -i "s/^NSTEPS .*/NSTEPS ${NSTEPS}/" "${SETTINGS}"
sed -i "s/^AUTO_EQUILIBRATION .*/AUTO_EQUILIBRATION 0/" "${SETTINGS}"
sed -i "s/^NBLOCKS_EQUIL .*/NBLOCKS_EQUIL ${NBLOCKS_EQUIL}/" "${SETTINGS}"
sed -i "s/^NSTEPS_EQUIL .*/NSTEPS_EQUIL ${NSTEPS_EQUIL}/" "${SETTINGS}"

sed -i "s/^NPART .*/NPART ${NPART}/" "${SETTINGS}"
for T in "${TEMPS[@]}"; do
  sed -i "s/^TEMP .*/TEMP ${T}/" "${SETTINGS}"
  for ALGO in "${ALGOS[@]}"; do
    sed -i "s/^ALGORITHM .*/ALGORITHM ${ALGO}/" "${IS_SETT}"

    DIR="${RES_DIR}/${ALGO}_T${T}"
    mkdir -p "${DIR}"

    sed -i "s/^H .*/H 0/" "${IS_SETT}"
    run
    rm "${OUTP_DIR}/magnetization.csv"
    mv "${OUTP_DIR}"/* "${DIR}"

    sed -i "s/^H .*/H ${H}/" "${IS_SETT}"
    run
    mv "${OUTP_DIR}/magnetization.csv" "${DIR}"
  done
done

