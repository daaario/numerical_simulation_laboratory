#!/bin/bash
set -euo pipefail

SETTINGS="input/sys_settings.dat"
LJ_SETT="input/lj_settings.dat"

OUTP_DIR="output"
RES_DIR="../nsl_exercises/sim_results/ex7"

rm -rf "${RES_DIR}"/*

run() {
  ./source/run
}

NPART=108
TEMPS=(0.8 1.1 1.2)
RHOS=(1.1 0.8 0.05)
PHASES=("SOLID" "LIQUID" "GAS")
RCUTS=(2.2 2.5 5.0)
NSTEPS=3600
NBLOCKS=16
NSTEPS_EQUIL=3600
NBLOCKS_EQUIL=12

sed -i "s/^MODEL .*/MODEL LJ/" "${SETTINGS}"
sed -i "s/^RESTARTED .*/RESTARTED 0/" "${SETTINGS}"
sed -i "s/^NPART .*/NPART ${NPART}/" "${SETTINGS}"
sed -i "s/^NBLOCKS .*/NBLOCKS ${NBLOCKS}/" "${SETTINGS}"
sed -i "s/^NSTEPS .*/NSTEPS ${NSTEPS}/" "${SETTINGS}"
sed -i "s/^AUTO_EQUILIBRATION .*/AUTO_EQUILIBRATION 1/" "${SETTINGS}"
sed -i "s/^NBLOCKS_EQUIL .*/NBLOCKS_EQUIL ${NBLOCKS_EQUIL}/" "${SETTINGS}"
sed -i "s/^NSTEPS_EQUIL .*/NSTEPS_EQUIL ${NSTEPS_EQUIL}/" "${SETTINGS}"

for i in "${!TEMPS[@]}"; do
  T="${TEMPS[i]}"
  RHO="${RHOS[i]}"
  RCUT="${RCUTS[i]}"
  PHASE="${PHASES[i]}"
  
  sed -i "s/^ALGORITHM .*/ALGORITHM VERLET/" "${LJ_SETT}"
  sed -i "s/^TEMP .*/TEMP ${T}/" "${SETTINGS}"
  sed -i "s/^RHO .*/RHO ${RHO}/" "${LJ_SETT}"
  sed -i "s/^RCUT .*/RCUT ${RCUT}/" "${LJ_SETT}"

  DIR="${RES_DIR}/VERLET_${PHASE}"
  mkdir -p "${DIR}"
  run
  T=$(cat "${OUTP_DIR}"/logfile.log | grep "^TEMP_EQUIL: " | awk '{printf $2}')
  mv "${OUTP_DIR}"/* "${DIR}"

  sed -i "s/^TEMP .*/TEMP ${T}/" "${SETTINGS}"
  sed -i "s/^ALGORITHM .*/ALGORITHM METROPOLIS/" "${LJ_SETT}"

  DIR="${RES_DIR}/METROPOLIS_${PHASE}"
  mkdir -p "${DIR}"
  run
  mv "${OUTP_DIR}"/* "${DIR}"
done


