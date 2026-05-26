#!/usr/bin/env bash
# Run mpboot placement (+ optional ratchet) on a GISAID test dataset.
#
# Usage:
#   ./run_gisaid.sh <test_size> <rep> [--ratchet] [--skip-build]
#
# <test_size> is one of: 9909, 19767, 29624, 42448
# <rep>       is one of: 0, 1, 2, 3, 4
#
# Examples:
#   ./run_gisaid.sh 9909 0                # placement + optimize, no ratchet (builds if needed)
#   ./run_gisaid.sh 9909 0 --ratchet      # + ratchet (iter=10, seed=42, runs=1)
#   ./run_gisaid.sh 9909 0 --skip-build   # skip the cmake/make step
#
# Outputs land next to the VCF (data/gisaid/test<N>/input/alignments/alignment_<rep>.vcf.{treefile,benchmark.json}).
# This script does NOT touch data/gisaid/test<N>/results/ — those hold published thesis numbers.

set -euo pipefail

USE_RATCHET=0
SKIP_BUILD=0
POSITIONAL=()

for arg in "$@"; do
  case "$arg" in
    --ratchet)    USE_RATCHET=1 ;;
    --skip-build) SKIP_BUILD=1 ;;
    --help|-h)
      sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'
      exit 0
      ;;
    -*)
      echo "Unknown flag: $arg" >&2
      exit 2
      ;;
    *)
      POSITIONAL+=("$arg")
      ;;
  esac
done

if [[ ${#POSITIONAL[@]} -ne 2 ]]; then
  echo "Usage: $0 <test_size> <rep> [--ratchet] [--skip-build]" >&2
  exit 2
fi

TEST_SIZE="${POSITIONAL[0]}"
REP="${POSITIONAL[1]}"

# Dataset-specific (pp_n = base/80%, pp_k = test/20%). Source of truth:
# thesis/benchmark_scripts/run_gisaid_sproptimizer.sh.
# Plain case stmt — macOS ships bash 3.2, no associative arrays.
case "$TEST_SIZE" in
  9909)  PP_N=7928;  PP_K=1981 ;;
  19767) PP_N=15814; PP_K=3953 ;;
  29624) PP_N=23700; PP_K=5924 ;;
  42448) PP_N=33959; PP_K=8489 ;;
  *)
    echo "Unknown GISAID test_size: $TEST_SIZE (expected one of: 9909, 19767, 29624, 42448)" >&2
    exit 2
    ;;
esac

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MPBOOT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$MPBOOT_ROOT/.." && pwd)"

BUILD_DIR="$MPBOOT_ROOT/build"
MPBOOT_BIN="$BUILD_DIR/mpboot"
DATA_DIR="$REPO_ROOT/data/gisaid/test${TEST_SIZE}"
VCF="$DATA_DIR/input/alignments/alignment_${REP}.vcf"
TREE="$DATA_DIR/input/trees/pruned_tree_${REP}.treefile"

if [[ ! -f "$VCF" ]]; then
  echo "VCF not found: $VCF" >&2
  exit 1
fi
if [[ ! -f "$TREE" ]]; then
  echo "Tree not found: $TREE" >&2
  exit 1
fi

if [[ "$SKIP_BUILD" -eq 0 ]]; then
  if [[ ! -d "$BUILD_DIR" ]]; then
    echo "+ mkdir $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
  fi
  if [[ ! -f "$BUILD_DIR/Makefile" && ! -f "$BUILD_DIR/build.ninja" ]]; then
    echo "+ cmake $MPBOOT_ROOT -DIQTREE_FLAGS=sse4 -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++"
    (cd "$BUILD_DIR" && cmake "$MPBOOT_ROOT" \
      -DIQTREE_FLAGS=sse4 \
      -DCMAKE_C_COMPILER=clang \
      -DCMAKE_CXX_COMPILER=clang++)
  fi
  # Rebuild if binary missing or any source is newer than the binary.
  NEEDS_BUILD=1
  if [[ -x "$MPBOOT_BIN" ]]; then
    NEWER=$(find "$MPBOOT_ROOT" \
      -path "$BUILD_DIR" -prune -o \
      -type f \( -name '*.cpp' -o -name '*.h' -o -name 'CMakeLists.txt' \) \
      -newer "$MPBOOT_BIN" -print -quit 2>/dev/null || true)
    [[ -z "$NEWER" ]] && NEEDS_BUILD=0
  fi
  if [[ "$NEEDS_BUILD" -eq 1 ]]; then
    echo "+ make -j4 (in $BUILD_DIR)"
    (cd "$BUILD_DIR" && make -j4)
  else
    echo "+ build up-to-date, skipping make"
  fi
fi

if [[ ! -x "$MPBOOT_BIN" ]]; then
  echo "mpboot binary missing after build: $MPBOOT_BIN" >&2
  exit 1
fi

CMD=(
  "$MPBOOT_BIN"
  -s "$VCF"
  -pp_tree "$TREE"
  -pp_on
  -pp_n "$PP_N"
  -pp_k "$PP_K"
  -pp_optimize
)

if [[ "$USE_RATCHET" -eq 1 ]]; then
  CMD+=(-pp_ratchet_iter 10 -pp_ratchet_seed 42 -pp_ratchet_runs 1)
fi

echo "+ ${CMD[*]}"
exec "${CMD[@]}"
