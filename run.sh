#!/bin/bash
set -e # exit immediately if a command exits with a non-zero status

echo "Creating virtual environment with Python..."
python -m venv venv

echo "activating virtual environment..."
souce venv/bin/activate

echo "Installing required Python packages..."
pip install --upgrade pip
pip install -r requirement.txt

echo "Cloning and installing tsconvert..."
git clone https://github.com/tskit-dev/tsconvert.git
cd tsconvert
pip install .
cd ..

echo "Cloning and building Relate..."
git clone https://github.com/MyersGroup/relate.git
cd relate/build

# Try to load CMake and GCC modules (if on an HPC or module system)
if command -v module &> /dev/null; then
    echo "Trying to load cmake and gcc modules..."
    if ! module load cmake/3.22.5 gcc/11.3.0; then
        echo "Module load failed. Attempting to install CMake locally..."
    fi
else
    echo "No module system detected. Checking if cmake is available..."
fi

# Check if cmake is installed
if ! command -v cmake &> /dev/null; then
    echo "CMake not found. Attempting to install..."

    # Try installing with apt (Debian/Ubuntu)
    if command -v apt &> /dev/null; then
        sudo apt update
        sudo apt install -y cmake
    # Try installing with Homebrew (macOS)
    elif command -v brew &> /dev/null; then
        brew install cmake
    else
        echo "Could not install CMake automatically. Please install it manually."
        exit 1
    fi
else
    echo "CMake is already available: $(cmake --version | head -n 1)"
fi

echo "Running cmake and make to build Relate..."
cmake..
make
cd ../..

echo "Setup complete. Running spacetrees via Snakemake with 1 thread..."
snakemake all -c1

echo "Done!"


