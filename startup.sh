#!/bin/bash
source venv/bin/activate
echo "Activated Virtual Env"
export XDG_CACHE_HOME=/scratch/m/mmosmond/raghavs/.cache
echo "Set Cache Location"
module load NiaEnv/2022a
module load gcc
echo "Loaded gcc"
gcc --version
