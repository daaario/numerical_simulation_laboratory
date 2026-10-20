#!/bin/bash
source ~/.george.env

# temperature_sweep.sh - Runs the simulation for an evenly-spaced sequence of temperatures, 
#                        saves results in 'output/T_*' for each temperature value.
# usage: ./temperature_sweep.sh <T_min> <T_max> <step>

set -euo pipefail

SETTINGS="input/sys_settings.dat"
RESULTS_DIR="output_tmp"
SETTINGS_BU="${SETTINGS}.backup.tmp"

run() {
  ./source/run
}

cleanup() {
  if [ -d ${RESULTS_DIR} ]; then
    if [ -d output ]; then 
      rm -rf output
    fi
    mv "${RESULTS_DIR}" output
  fi

  if [ -f ${SETTINGS_BU} ]; then
    cp "${SETTINGS_BU}" "${SETTINGS}"
    rm "${SETTINGS_BU}"
  fi
}

trap cleanup EXIT INT TERM

if [ $# -eq 0 ]; then
  TEMPS=(0.1 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0)
  curl \
    --connect-timeout 5 \
    -H "Priority: 5" \
    -H "Tags: gear" \
    -H "Title: Simulation started." \
    -d "Running temperature sweep over 0.1, 0.5, 1.0, 1.5, ..., 5.0." \
    "$NTFY_URL" 1>/dev/null 2>&1
elif [ $# -eq 3 ]; then
  TEMPS=($(seq $1 $3 $2))
  curl \
    --connect-timeout 5 \
    -H "Priority: 5" \
    -H "Tags: gear" \
    -H "Title: Simulation started." \
    -d "Running temperature sweep from T=$1 to T=$2, with step dT=$3." \
    "$NTFY_URL" 1>/dev/null 2>&1
else
  echo "Usage:  ./temperature_sweep.sh <T_min> <T_max> <step>
  If no arguments are provided, the default temperatures are 0.1, 0.5, 1.0, 1.5, ..., 5.0."
  exit 1
fi

make 1>/dev/null 2>&1
mkdir -p "$RESULTS_DIR"
cp "$SETTINGS" "${SETTINGS_BU}" 

for T in "${TEMPS[@]}"; do
  echo "  ============== Running sim at T=${T} =============="
  curl --connect-timeout 5 -d "Computing for T=${T}" "${NTFY_URL}" 1>/dev/null 2>&1

  # from the template, edit 'output/sys.settings' according to each temperature
  sed "s/^TEMP .*/TEMP $T/" "${SETTINGS_BU}" > "$SETTINGS"
  sed -i "s/^RESTARTED .*/RESTARTED 0/" "$SETTINGS"

  T_DIR="$RESULTS_DIR/T_${T}"
  rm -rf output
  mkdir -p "$T_DIR"
  ln -s "$(realpath "${T_DIR}")" output

  run
done

curl \
  --connect-timeout 5 \
  -H "Priority: 5" \
  -H "Tags: coffee" \
  -H "Title: Simulation finished!" \
  -d "Results are avaiable in 'output' directory." \
  "$NTFY_URL" 1>/dev/null 2>&1
