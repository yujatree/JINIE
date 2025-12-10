#!/bin/bash

# Colors
RED="\033[0;31m"
GREEN="\033[0;32m"
BLUE="\033[0;34m"
NC="\033[0m"   # No Color

echo
echo "===================================================="
echo
echo "        █████ █████ ██████   █████ █████ ██████████"
echo "       ░░███ ░░███ ░░██████ ░░███ ░░███ ░░███░░░░░█"
echo "        ░███  ░███  ░███░███ ░███  ░███  ░███  █ ░ "
echo "        ░███  ░███  ░███░░███░███  ░███  ░██████   "
echo "        ░███  ░███  ░███ ░░██████  ░███  ░███░░█   "
echo "  ███   ░███  ░███  ░███  ░░█████  ░███  ░███ ░   █"
echo " ░░████████   █████ █████  ░░█████ █████ ██████████"
echo "  ░░░░░░░░   ░░░░░ ░░░░░    ░░░░░ ░░░░░ ░░░░░░░░░░ "
echo " --------------------------------------------------"
echo "          ▗ ▌     ▖ ▄▖▄▖        ▜     ▘    ▗     ▜ " 
echo "      ▄▖  ▜▘▛▌█▌  ▌ ▙▌▚   ▀▌▛▌▀▌▐ ▌▌▛▘▌▛▘  ▜▘▛▌▛▌▐ "
echo "          ▐▖▌▌▙▖  ▙▖▌ ▄▌  █▌▌▌█▌▐▖▙▌▄▌▌▄▌  ▐▖▙▌▙▌▐▖"
echo "                                  ▄▌               "                        
echo "===================================================="
echo
echo
echo " [1] Mean Squared Displacement"
echo " [2] Non-Gaussian Parameter"
echo " [3] Ionic Conductivity"
echo " [4] Radial Distribution Function"
echo " [5] Self-part van Hove Correlation Function"
echo " [6] Hop Function"
echo " [7] DCD to Lammps trajectory"
echo " [8] EXIT"
echo
read -p " [Enter the property index] : " idx

case $idx in
    1) prop="MSD" ;;
    2) prop="NGP" ;;
    3) prop="ICD" ;;
    4) prop="RDF" ;;
    5) prop="GRT" ;;
    6) prop="HOP" ;;
    7) prop="D2L" ;;
    8) echo ; echo " Program ends ..."; echo 
       echo "====================================================" ; echo
       exit 0 
       ;;
    *) echo " ERROR: Enter the right index." ; echo
       echo "====================================================" ; echo
       exit 1 
       ;;
esac

# Execute file
execfile="${HOME}/CODES/JINIE/bin/${prop}"
if [[ ! -x "$execfile" ]]; then
    echo " ERROR: There is no exec. file or no permission: "
    echo "             $execfile"
    echo ; echo "===================================================="; echo
    exit 1
fi

# Input trajectories
echo 
read -p " Trajectories for analysis: " -a inputs   # 배열로 입력받음

if [[ ${#inputs[@]} -eq 0 ]]; then
    echo " ERROR: You have to select at least one trajectory file." ; echo
    echo "====================================================" ; echo
    exit 1
fi

# Inspectation of extensions
first_ext="${inputs[0]##*.}" # Set standard extension

if [[ "$first_ext" != "dcd" && "$first_ext" != "lammpstraj" ]]; then
    echo " ERROR: Only supports .dcd or ./lammpstraj." ; echo
    echo "====================================================" ; echo
    exit 1
fi

for f in "${inputs[@]}"; do
    if [[ ! -f "$f" ]]; then
        echo " ERROR: There is no file.: $f" ; echo
        echo "====================================================" ; echo
        exit 1
    fi

    ext="${f##*.}"
    if [[ "$ext" != "$first_ext" ]]; then
        echo " ERROR: All of the input files should have the same extension."
        echo " - Standard: .$first_ext"
        echo " - One of the files: $f (.$ext)" ; echo
        echo "====================================================" ; echo
        exit 1
    fi
done

# log directory
mkdir -p ./log
logfile="./log/${prop}.log"

# execute data in log file
start_time=$(date "+%Y-%m-%d %H:%M:%S")
echo "====== $prop Calculation start ======" > "$logfile"
echo "Start   : $start_time" >> "$logfile"
echo "Inputs  : ${inputs[*]}" >> "$logfile"
echo "===================================" >> "$logfile"

# EXECUTE
echo 
nohup "$execfile" "${inputs[@]}" >> "$logfile" 2>&1 &

echo " Program : $prop"
echo " PID     : $pid"
echo " Log     : $logfile"
echo 
echo "===================================================="
