#!/bin/bash
#
# MPBoot Optimization Benchmark Script (SPR + Ratchet + TBR)
#
# Tests each case with multiple SPR radius values: max(32, x% of test_size)
# for x = 5, 7, 9, 11 (configurable via --radius_pcts)
#
# Usage:
#   ./scripts/run_mpboot_benchmark.sh                          # Run all tests
#   ./scripts/run_mpboot_benchmark.sh 200                      # Run all indices for test200
#   ./scripts/run_mpboot_benchmark.sh 200 1                    # Run only test200/1
#   ./scripts/run_mpboot_benchmark.sh --passes 3               # Run all with 3 SPR passes
#   ./scripts/run_mpboot_benchmark.sh 200 1 --passes 3         # Specific test with 3 passes
#   ./scripts/run_mpboot_benchmark.sh --pp_k 300 --pp_n 300    # Override scaled defaults
#   ./scripts/run_mpboot_benchmark.sh --radius_pcts "3,5,7,9"  # Custom radius percentages
#
# pp_k/pp_n default to test size (overridable with --pp_k/--pp_n)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SOURCE_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ROOT_DIR="$(cd "$SOURCE_DIR/.." && pwd)"

# Support override via env vars (for Docker)
BUILD_DIR="${MPBOOT_BUILD_DIR:-$ROOT_DIR/build}"
DATA_DIR="${MPBOOT_DATA_DIR:-$ROOT_DIR/placement_genbank_data}"
MPBOOT="$BUILD_DIR/mpboot"

# Defaults
TARGET_SIZE=""
TARGET_INDEX=""
SPR_PASSES=1
PP_K_OVERRIDE=""
PP_N_OVERRIDE=""
RADIUS_PCTS="5,7,9,11"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --passes)      SPR_PASSES="$2"; shift 2 ;;
        --pp_k)        PP_K_OVERRIDE="$2"; shift 2 ;;
        --pp_n)        PP_N_OVERRIDE="$2"; shift 2 ;;
        --radius_pcts) RADIUS_PCTS="$2"; shift 2 ;;
        *)
            if [[ -z "$TARGET_SIZE" ]]; then
                TARGET_SIZE="$1"
            elif [[ -z "$TARGET_INDEX" ]]; then
                TARGET_INDEX="$1"
            fi
            shift ;;
    esac
done

# Convert comma-separated percentages to array
IFS=',' read -ra PCT_ARRAY <<< "$RADIUS_PCTS"

# Compute radius values for a given test size: max(32, size * pct / 100)
compute_radii() {
    local size=$1
    local radii=()
    local seen=()
    for pct in "${PCT_ARRAY[@]}"; do
        local r=$(( size * pct / 100 ))
        if (( r < 32 )); then
            r=32
        fi
        # Deduplicate (small tests may produce the same radius for multiple pcts)
        local dup=0
        for s in "${seen[@]+"${seen[@]}"}"; do
            if (( s == r )); then dup=1; break; fi
        done
        if (( dup == 0 )); then
            radii+=("$r")
            seen+=("$r")
        fi
    done
    echo "${radii[@]}"
}

# pp_k and pp_n default to the test size (matching run_test.sh behavior).
get_scaled_pp() {
    echo "$1"
}

# Verify build exists
if [[ ! -x "$MPBOOT" ]]; then
    echo "ERROR: mpboot binary not found at $MPBOOT"
    echo "Build first: cd $BUILD_DIR && cmake ../mpboot -DCMAKE_BUILD_TYPE=Release && make -j\$(nproc)"
    exit 1
fi

# CSV output (override with MPBOOT_OUTPUT_DIR env var)
OUTPUT_DIR="${MPBOOT_OUTPUT_DIR:-$SOURCE_DIR}"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
CSV_FILE="$OUTPUT_DIR/mpboot-optimize-result-${TIMESTAMP}.csv"
echo "test_size,index,pp_k,pp_n,spr_passes,spr_max_radius,radius_pct,placement_score,placement_time_s,spr_score,spr_time_s,fitch_verify_score,delta,total_time_s,status" > "$CSV_FILE"

# Discover test sizes
if [[ -n "$TARGET_SIZE" ]]; then
    SIZES=("$TARGET_SIZE")
else
    SIZES=()
    for dir in "$DATA_DIR"/test*; do
        size=$(basename "$dir" | sed 's/test//')
        SIZES+=("$size")
    done
    # Sort numerically
    IFS=$'\n' SIZES=($(sort -n <<<"${SIZES[*]}")); unset IFS
fi

TOTAL_PASS=0
TOTAL_FAIL=0
TOTAL_SKIP=0

for SIZE in "${SIZES[@]}"; do
    TEST_DIR="$DATA_DIR/test${SIZE}"
    if [[ ! -d "$TEST_DIR" ]]; then
        echo "WARNING: $TEST_DIR not found, skipping"
        continue
    fi

    # Discover indices
    if [[ -n "$TARGET_INDEX" ]]; then
        INDICES=("$TARGET_INDEX")
    else
        INDICES=()
        for dir in "$TEST_DIR"/*/; do
            idx=$(basename "$dir")
            [[ "$idx" =~ ^[0-9]+$ ]] && INDICES+=("$idx")
        done
        IFS=$'\n' INDICES=($(sort -n <<<"${INDICES[*]}")); unset IFS
    fi

    # Determine pp_k/pp_n for this test size (CLI override takes precedence)
    if [[ -n "$PP_K_OVERRIDE" ]]; then
        PP_K="$PP_K_OVERRIDE"
    else
        PP_K=$(get_scaled_pp "$SIZE")
    fi
    if [[ -n "$PP_N_OVERRIDE" ]]; then
        PP_N="$PP_N_OVERRIDE"
    else
        PP_N=$(get_scaled_pp "$SIZE")
    fi

    # Compute radius values for this test size
    RADII=($(compute_radii "$SIZE"))
    echo "=== test${SIZE}: radii = ${RADII[*]} (from pcts: ${RADIUS_PCTS}) ==="

    for IDX in "${INDICES[@]}"; do
        CASE_DIR="$TEST_DIR/$IDX"
        VCF="$CASE_DIR/added${SIZE}.vcf"
        TREE="$CASE_DIR/origin${SIZE}.fasta.treefile"

        if [[ ! -f "$VCF" ]]; then
            echo "SKIP test${SIZE}/${IDX}: VCF not found ($VCF)"
            echo "${SIZE},${IDX},${PP_K},${PP_N},${SPR_PASSES},,,,,,,,skip_no_vcf" >> "$CSV_FILE"
            ((TOTAL_SKIP++))
            continue
        fi
        if [[ ! -f "$TREE" ]]; then
            echo "SKIP test${SIZE}/${IDX}: tree not found ($TREE)"
            echo "${SIZE},${IDX},${PP_K},${PP_N},${SPR_PASSES},,,,,,,,skip_no_tree" >> "$CSV_FILE"
            ((TOTAL_SKIP++))
            continue
        fi

        for RADIUS_IDX in "${!RADII[@]}"; do
            RADIUS="${RADII[$RADIUS_IDX]}"
            # Figure out which percentage produced this radius
            RADIUS_PCT="${PCT_ARRAY[$RADIUS_IDX]}"

            echo "--- Running test${SIZE}/${IDX} (radius=${RADIUS}, ${RADIUS_PCT}%, pp_k=${PP_K}, pp_n=${PP_N}, passes=${SPR_PASSES}) ---"

            TMPLOG=$(mktemp)
            START_TS=$(date +%s)

            set +e
            "$MPBOOT" \
                -s "$VCF" \
                -pp_tree "$TREE" \
                -pp_on \
                -pp_k "$PP_K" \
                -pp_n "$PP_N" \
                -mpboot_optimize \
                -spr_max_passes "$SPR_PASSES" \
                -spr_max_radius "$RADIUS" \
                > "$TMPLOG" 2>&1
            EXIT_CODE=$?
            set -e

            END_TS=$(date +%s)
            ELAPSED=$((END_TS - START_TS))

            # Check for errors in output (mpboot sometimes exits 0 on error)
            if grep -q "ERROR:" "$TMPLOG" 2>/dev/null; then
                ERROR_MSG=$(grep "ERROR:" "$TMPLOG" | head -1 | sed 's/.*ERROR: //')
                echo "  ERROR: ${ERROR_MSG} (${ELAPSED}s)"
                echo "${SIZE},${IDX},${PP_K},${PP_N},${SPR_PASSES},${RADIUS},${RADIUS_PCT},,,,,,${ELAPSED},fail_error" >> "$CSV_FILE"
                ((TOTAL_FAIL++))
            elif [[ $EXIT_CODE -ne 0 ]]; then
                echo "  CRASHED (exit code ${EXIT_CODE}) ${ELAPSED}s"
                echo "${SIZE},${IDX},${PP_K},${PP_N},${SPR_PASSES},${RADIUS},${RADIUS_PCT},,,,,,${ELAPSED},fail_crash" >> "$CSV_FILE"
                ((TOTAL_FAIL++))
            else
                # Parse scores
                PLACEMENT_SCORE=$(grep "Placement parsimony score" "$TMPLOG" | grep -oE '[0-9]+' | tail -1 || echo "")
                SPR_SCORE=$(grep "Final parsimony score after SPR" "$TMPLOG" | grep -oE '[0-9]+' | tail -1 || echo "")
                FITCH_SCORE=$(grep "Final parsimony score computed by fitch" "$TMPLOG" | grep -oE '[0-9]+' | tail -1 || echo "")

                # Parse times (e.g. "Time: 0.020 seconds", "SPR optimization time: 3.782 seconds")
                PLACEMENT_TIME=$(grep "Finished placement core" -A1 "$TMPLOG" | grep "Time:" | grep -oE '[0-9]+\.[0-9]+' | head -1 || echo "")
                SPR_TIME=$(grep "SPR optimization time:" "$TMPLOG" | grep -oE '[0-9]+\.[0-9]+' | head -1 || echo "")

                if [[ -n "$PLACEMENT_SCORE" && -n "$SPR_SCORE" ]]; then
                    DELTA=$((PLACEMENT_SCORE - SPR_SCORE))
                else
                    DELTA=""
                fi

                # Verify fitch matches SPR
                if [[ -n "$SPR_SCORE" && -n "$FITCH_SCORE" && "$SPR_SCORE" == "$FITCH_SCORE" ]]; then
                    STATUS="pass"
                    ((TOTAL_PASS++))
                elif [[ -n "$SPR_SCORE" && -n "$FITCH_SCORE" ]]; then
                    STATUS="fail_mismatch"
                    ((TOTAL_FAIL++))
                else
                    STATUS="fail_parse"
                    ((TOTAL_FAIL++))
                fi

                echo "  Placement: ${PLACEMENT_SCORE:-?} (${PLACEMENT_TIME:-?}s) -> SPR: ${SPR_SCORE:-?} (${SPR_TIME:-?}s) radius=${RADIUS} delta=${DELTA:-?} fitch=${FITCH_SCORE:-?} [${STATUS}] total=${ELAPSED}s"
                echo "${SIZE},${IDX},${PP_K},${PP_N},${SPR_PASSES},${RADIUS},${RADIUS_PCT},${PLACEMENT_SCORE},${PLACEMENT_TIME},${SPR_SCORE},${SPR_TIME},${FITCH_SCORE},${DELTA},${ELAPSED},${STATUS}" >> "$CSV_FILE"
            fi

            rm -f "$TMPLOG"
        done
    done
done

echo ""
echo "========================================="
echo "Results: ${TOTAL_PASS} passed, ${TOTAL_FAIL} failed, ${TOTAL_SKIP} skipped"
echo "CSV: $CSV_FILE"
echo "========================================="
