#!/bin/bash

# Define variables with default values
SCRIPT_PATH="/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/glycoclust/scripts/glycoclust.R"
CONFIG_PATH=""
SEED=123456

# Print usage information
usage() {
    echo "Usage: $0 -c <config_path> [-s <seed>]"
    echo "  -c  Path to the configuration file (required)"
    echo "  -s  Random seed (optional, default: 123456)"
    exit 1
}

# Parse command-line arguments
while getopts "c:s:" opt; do
    case $opt in
        c) CONFIG_PATH="$OPTARG" ;;
        s) SEED="$OPTARG" ;;
        *) usage ;;
    esac
done

# Check if config path is provided
if [ -z "$CONFIG_PATH" ]; then
    echo "Error: Configuration path (-c) is required."
    usage
fi

# Run the glycoclust script
Rscript "$SCRIPT_PATH" -c "$CONFIG_PATH" -s "$SEED"
