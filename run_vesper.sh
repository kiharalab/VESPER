#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="${SCRIPT_DIR}/.venv"

if [[ ! -d "$VENV_DIR" ]]; then
    echo "=== Creating virtual environment ==="
    python3 -m venv "$VENV_DIR"
    "$VENV_DIR/bin/pip" install --upgrade pip
    "$VENV_DIR/bin/pip" install numpy scipy mrcfile
fi

source "$VENV_DIR/bin/activate"

usage() {
    cat <<EOF
Usage: $0 -a MAP1.mrc -b MAP2.mrc [VESPER options...]

Unifies both input maps (fixing axis order and origin) before running VESPER.

Required:
  -a MAP1.mrc   First (larger) EM map
  -b MAP2.mrc   Second (smaller) EM map

All remaining arguments are passed directly to VESPER.

Example:
  $0 -a emd_8724.map -b emd_8409.map -t 0.04 -T 0.048 -s 7 -A 30 -c 5 -S
EOF
    exit 1
}

# Parse -a and -b from arguments, collect the rest for VESPER
MAP1=""
MAP2=""
VESPER_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -a)
            MAP1="$2"
            shift 2
            ;;
        -b)
            MAP2="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            VESPER_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -z "$MAP1" || -z "$MAP2" ]]; then
    echo "Error: Both -a and -b map files are required."
    usage
fi

if [[ ! -f "$MAP1" ]]; then
    echo "Error: Map file not found: $MAP1"
    exit 1
fi

if [[ ! -f "$MAP2" ]]; then
    echo "Error: Map file not found: $MAP2"
    exit 1
fi

UNIFY_PY="${SCRIPT_DIR}/unify.py"
VESPER_BIN="${SCRIPT_DIR}/VESPER_code/VESPER"

if [[ ! -f "$VESPER_BIN" ]]; then
    echo "=== VESPER binary not found, compiling from source ==="
    make -C "${SCRIPT_DIR}/VESPER_code" || {
        echo "Error: Compilation failed. Make sure gcc and fftw3 are installed."
        exit 1
    }
    cp "${SCRIPT_DIR}/VESPER_code/VESPER" "${SCRIPT_DIR}/VESPER"
fi

# Generate unified map paths next to the originals
make_unified_path() {
    local dir
    local base
    local ext
    dir="$(dirname "$1")"
    base="$(basename "$1")"
    ext="${base##*.}"
    base="${base%.*}"
    echo "${dir}/${base}_unified.${ext}"
}

UNIFIED_MAP1="$(make_unified_path "$MAP1")"
UNIFIED_MAP2="$(make_unified_path "$MAP2")"

echo "=== Unifying Map 1: $MAP1 -> $UNIFIED_MAP1 ==="
python3 "$UNIFY_PY" -i "$MAP1" -o "$UNIFIED_MAP1"

echo ""
echo "=== Unifying Map 2: $MAP2 -> $UNIFIED_MAP2 ==="
python3 "$UNIFY_PY" -i "$MAP2" -o "$UNIFIED_MAP2"

echo ""
echo "=== Running VESPER ==="
echo "$VESPER_BIN" -a "$UNIFIED_MAP1" -b "$UNIFIED_MAP2" "${VESPER_ARGS[@]}"
"$VESPER_BIN" -a "$UNIFIED_MAP1" -b "$UNIFIED_MAP2" "${VESPER_ARGS[@]}"
