#!/usr/bin/env bash

# Stop execution if any command fails
set -e

# Resolve the absolute path of this script, handling symlinks safely
TARGET_FILE=$0
while [ -L "$TARGET_FILE" ]; do
    TARGET_FILE=$(readlink "$TARGET_FILE")
done
REPO_DIR=$(cd "$(dirname "$TARGET_FILE")" && pwd)
MAIN_NF="${REPO_DIR}/main.nf"
CONFIG_FILE="${REPO_DIR}/nextflow.config"

# Sanity checks
if ! command -v nextflow &> /dev/null; then
    echo "Error: 'nextflow' is not installed or not found in your PATH." >&2
    exit 1
fi

if [ ! -f "$MAIN_NF" ]; then
    echo "Error: Cannot find main.nf at $MAIN_NF" >&2
    exit 1
fi

# Initialize variables for custom flags
FASTQ=""
FASTQ_1=""
FASTQ_2=""
PREFIX=""
SAMPLE_SHEET=""
NF_ARGS=()

# Helper function to convert paths to absolute paths and validate existence
get_abs_path() {
    if [[ -n "$1" ]]; then
        if [[ ! -e "$1" ]]; then
            echo "Error: File '$1' does not exist." >&2
            exit 1
        fi
        echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
    fi
}

# Parse custom wrapper flags, leave everything else for Nextflow
while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)
            if [[ -f "$CONFIG_FILE" ]]; then
                # Pulls out the text inside the single or double quotes of the version line
                VERSION=$(sed -nE "s/^[[:space:]]*version[[:space:]]*=[[:space:]]*['\"]([^'\"]+)['\"].*/\1/p" "$CONFIG_FILE")
                if [[ -n "$VERSION" ]]; then
                    echo "donut_falls version $VERSION"
                else
                    echo "Error: Could not parse version string from nextflow.config" >&2
                    exit 1
                fi
            else
                echo "Error: Cannot find nextflow.config at $CONFIG_FILE" >&2
                exit 1
            fi
            exit 0
            ;;
        --fastq)
            FASTQ="$2"
            shift 2
            ;;
        --fastq_1)
            FASTQ_1="$2"
            shift 2
            ;;
        --fastq_2)
            FASTQ_2="$2"
            shift 2
            ;;
        --prefix)
            PREFIX="$2"
            shift 2
            ;;
        --sample_sheet)
            SAMPLE_SHEET="$2"
            shift 2
            ;;
        *)
            NF_ARGS+=("$1")
            shift
            ;;
    esac
done

# --- Routing logic ---

# Conflict Check: Make sure they didn't pass both entrypoints
if [[ -n "$SAMPLE_SHEET" && -n "$FASTQ" ]]; then
    echo "Error: Cannot combine direct file inputs (--fastq) with an explicit --sample_sheet." >&2
    exit 1
fi

# Mode A: User supplied an explicit sample sheet
if [[ -n "$SAMPLE_SHEET" ]]; then
    ABS_SHEET=$(get_abs_path "$SAMPLE_SHEET")
    NF_ARGS+=("--sample_sheet" "$ABS_SHEET")
    
    if [[ -n "$PREFIX" ]]; then
        echo "Warning: The --prefix flag is ignored when running with a --sample_sheet." >&2
    fi

# Mode B: User supplied direct files (Wrapper generates temporary sample sheet)
elif [[ -n "$FASTQ" ]]; then
    if [[ -n "$PREFIX" ]]; then
        SAMPLE_ID="$PREFIX"
    else
        BASE_NAME=$(basename "$FASTQ")
        SAMPLE_ID=$(echo "$BASE_NAME" | cut -d. -f1)
    fi

    ABS_FASTQ=$(get_abs_path "$FASTQ")
    ABS_FASTQ_1=$(get_abs_path "$FASTQ_1")
    ABS_FASTQ_2=$(get_abs_path "$FASTQ_2")

    TEMP_SHEET=$(mktemp /tmp/donut_falls_sheet.XXXXXX.csv)
    trap 'rm -f "$TEMP_SHEET"' EXIT

    echo "sample,fastq,fastq_1,fastq_2" > "$TEMP_SHEET"
    echo "${SAMPLE_ID},${ABS_FASTQ},${ABS_FASTQ_1},${ABS_FASTQ_2}" >> "$TEMP_SHEET"

    NF_ARGS+=("--sample_sheet" "$TEMP_SHEET")

# Catch standalone flag errors
elif [[ -n "$PREFIX" ]]; then
    echo "Error: The --prefix flag can only be used alongside --fastq file inputs." >&2
    exit 1
fi

# Transparently execute Nextflow with the cleaned array
exec nextflow run "$MAIN_NF" "${NF_ARGS[@]}"