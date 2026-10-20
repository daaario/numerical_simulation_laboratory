#!/bin/bash

# ./plot <target_directory>
# Detects whether <target_directory> holds Ising or Lennard-Jones results by
# reading the first word of logfile.log in one of its T_* subdirectories, then
# executes the matching pair of plotting scripts in the background.
# If no target directory is specified, 'output' is used.

if [ $# -eq 0 ]; then
  DIR="output"
elif [ $# -eq 1 ]; then
  DIR=$(echo ${1} | sed 's/\/*$//')
else
  echo "
./plot <target_directory>
Detects whether <target_directory> holds Ising or Lennard-Jones results, then
executes the matching pair of plotting scripts in the background.
If no target directory is specified, 'output' is used.
"
  exit 1
fi

# grab logfile.log from any one T_* subdirectory -- system type is the same for all of them
LOGFILE=$(find "${DIR}" -maxdepth 2 -path "${DIR}/T_*/logfile.log" -print -quit)

if [ -z "${LOGFILE}" ]; then
  echo "Error: no T_*/logfile.log found under '${DIR}'." >&2
  exit 1
fi

SYSTEM=$(head -n 1 "${LOGFILE}" | awk '{print $1}')

case "${SYSTEM}" in
  Ising)
    ./source/plot_ising.py "${DIR}" 2>/dev/null & disown
#    ./source/plot_blocks_deviation.py "${DIR}" 2>/dev/null & disown
    ;;
  Lennard-Jones)
    ./source/plot_lj.py "${DIR}" 2>/dev/null & disown
    ./source/plot_pofv.py "${DIR}" 2>/dev/null & disown
    ./source/plot_gofr.py "${DIR}" 2>/dev/null & disown
    ;;
  *)
    echo "Error: unrecognized system '${SYSTEM}' in ${LOGFILE} (expected 'Ising' or 'Lennard-Jones')." >&2
    exit 1
    ;;
esac
