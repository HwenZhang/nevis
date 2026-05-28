#!/bin/bash
# ============================================================
#  NEVIS HPC-style Job Runner
#  - File-based job state (portable, debuggable: see .jobs/)
#  - Auto retry on failure
#  - Timestamped per-job logs
#  - Spinup to drainage dependency chain
#  - Clean signal handling
#  - Works on macOS bash 3.x and Linux bash 4+
# ============================================================
set -uo pipefail

# --- Configuration (override via environment) ---
MAX_PARALLEL=${MAX_PARALLEL:-8}       # concurrent MATLAB jobs
MAX_RETRIES=${MAX_RETRIES:-2}         # retry failed jobs up to N times
SKIP_SPINUP=${SKIP_SPINUP:-false}      # skip spinup, run drainage directly
POLL_INTERVAL=${POLL_INTERVAL:-5}     # seconds between reap cycles
JOB_TIMEOUT=${JOB_TIMEOUT:-0}         # per-job wall-clock limit (0 = none)
CASE_DIR=${CASE_DIR:-cases}           # directory containing case configs
CASE_WORKFLOW=${CASE_WORKFLOW:-idealized} # idealized | regional

# --- Directories ---
LOG_DIR="logs"
STATE_DIR=".jobs"
DONE_DIR="parameter_sweep"
QUEUE_FILE="$STATE_DIR/queue"

# --- Counters ---
total_jobs=0
completed_jobs=0
failed_jobs=0

# ============================================================
#  Logging
# ============================================================
_ts() { date '+%Y-%m-%d %H:%M:%S'; }
log()      { printf '[%s]      %s\n' "$(_ts)" "$*"; }
log_ok()   { printf '[%s]  OK  %s\n' "$(_ts)" "$*"; }
log_warn() { printf '[%s] WARN %s\n' "$(_ts)" "$*"; }
log_err()  { printf '[%s] FAIL %s\n' "$(_ts)" "$*" >&2; }

# ============================================================
#  File-based state helpers
#  State dir layout:
#    .jobs/<name>.state    PENDING | RUNNING | DONE | FAILED | KILLED
#    .jobs/<name>.pid      PID of running MATLAB process
#    .jobs/<name>.retries  number of retries consumed so far
#    .jobs/<name>.start    epoch second when last launched
# ============================================================
state_set()   { echo "$2" > "$STATE_DIR/${1}.state"; }
state_get()   { cat "$STATE_DIR/${1}.state" 2>/dev/null || echo "NONE"; }
retries_get() { cat "$STATE_DIR/${1}.retries" 2>/dev/null || echo "0"; }
retries_set() { echo "$2" > "$STATE_DIR/${1}.retries"; }
pid_set()     { echo "$2" > "$STATE_DIR/${1}.pid"; }
pid_get()     { cat "$STATE_DIR/${1}.pid" 2>/dev/null || echo ""; }
start_set()   { echo "$2" > "$STATE_DIR/${1}.start"; }
start_get()   { cat "$STATE_DIR/${1}.start" 2>/dev/null || echo "0"; }

# ============================================================
#  Queue (one job name per line in a file)
# ============================================================
queue_push() { echo "$1" >> "$QUEUE_FILE"; }

queue_pop() {
    # Sets REPLY to next job name; returns 1 if empty.
    if [ ! -s "$QUEUE_FILE" ]; then
        REPLY=""
        return 1
    fi
    REPLY=$(head -1 "$QUEUE_FILE")
    # portable: temp-file approach works on both macOS and Linux
    tail -n +2 "$QUEUE_FILE" > "$QUEUE_FILE.tmp" && mv "$QUEUE_FILE.tmp" "$QUEUE_FILE"
    return 0
}

queue_size() {
    if [ -s "$QUEUE_FILE" ]; then
        wc -l < "$QUEUE_FILE" | tr -d ' '
    else
        echo 0
    fi
}

case_matches_workflow() {
    local name="$1"
    case "$CASE_WORKFLOW" in
        regional)
            case "$name" in
                n2d_region*|n2d_regional*) return 0 ;;
                *) return 1 ;;
            esac
            ;;
        idealized)
            case "$name" in
                n2d_region*|n2d_regional*) return 1 ;;
                n1d*|n2d*) return 0 ;;
                *) return 1 ;;
            esac
            ;;
        *)
            log_err "Unknown CASE_WORKFLOW=$CASE_WORKFLOW (expected idealized or regional)"
            exit 1
            ;;
    esac
}

# ============================================================
#  Running-job helpers (derived from state dir)
# ============================================================
get_running_names() {
    local names=""
    for f in "$STATE_DIR"/*.state; do
        [ -f "$f" ] || continue
        if [ "$(cat "$f")" = "RUNNING" ]; then
            local n
            n=$(basename "${f%.state}")
            names="$names $n"
        fi
    done
    echo "$names"
}

running_count() {
    local c=0
    for _ in $(get_running_names); do
        c=$((c + 1))
    done
    echo "$c"
}

# ============================================================
#  Cleanup (signal handler)
# ============================================================
cleanup() {
    echo ""
    log_warn "Signal received - stopping all running jobs..."
    for name in $(get_running_names); do
        local pid
        pid=$(pid_get "$name")
        if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null || true
            sleep 1
            kill -9 "$pid" 2>/dev/null || true
            log "  Killed $name (PID $pid)"
            state_set "$name" "KILLED"
        fi
    done
    # Safety net for any orphaned MATLAB processes
    pkill -f 'matlab.*-batch' 2>/dev/null || true
    log "Cleanup complete."
    print_summary
    exit 130
}
trap cleanup INT TERM

# ============================================================
#  Launch a single job
# ============================================================
launch_job() {
    local name="$1"
    local attempt
    attempt=$(retries_get "$name")
    local logfile="${LOG_DIR}/${name}.log"

    if [ "$attempt" -gt 0 ]; then
        printf '\n\n===== RETRY %d/%d  %s =====\n\n' "$attempt" "$MAX_RETRIES" "$(_ts)" >> "$logfile"
        log "Starting $name  [retry $attempt/$MAX_RETRIES]"
    else
        log "Starting $name"
    fi

    state_set "$name" "RUNNING"
    start_set "$name" "$(date +%s)"

    # Launch MATLAB; stdout+stderr -> per-job log
    if [ "$CASE_WORKFLOW" = "regional" ]; then
        matlab -batch "n2d_regional_template('$name')" >> "$logfile" 2>&1 &
    else
        matlab -batch "n2d_idealized_template('$name')" >> "$logfile" 2>&1 &
    fi
    local pid=$!
    pid_set "$name" "$pid"

    log "  PID=$pid  log=$logfile"
}

# ============================================================
#  Reap finished / timed-out jobs
# ============================================================
reap_jobs() {
    local now
    now=$(date +%s)

    for name in $(get_running_names); do
        local pid
        pid=$(pid_get "$name")
        [ -z "$pid" ] && continue

        # --- timeout check ---
        if [ "$JOB_TIMEOUT" -gt 0 ] && kill -0 "$pid" 2>/dev/null; then
            local started elapsed
            started=$(start_get "$name")
            elapsed=$((now - started))
            if [ "$elapsed" -gt "$JOB_TIMEOUT" ]; then
            log_warn "$name timed out (${elapsed}s > ${JOB_TIMEOUT}s) - killing"
                kill "$pid" 2>/dev/null || true
                sleep 1
                kill -9 "$pid" 2>/dev/null || true
                state_set "$name" "FAILED"
                completed_jobs=$((completed_jobs + 1))
                failed_jobs=$((failed_jobs + 1))
                log_err "$name  (timeout after ${elapsed}s)"
                continue
            fi
        fi

        # --- still running? ---
        if kill -0 "$pid" 2>/dev/null; then
            continue
        fi

        # --- finished: collect exit code ---
        wait "$pid" 2>/dev/null
        local exit_code=$?
        local started elapsed
        started=$(start_get "$name")
        elapsed=$((now - started))

        if [ "$exit_code" -eq 0 ]; then
            # ---- success ----
            state_set "$name" "DONE"
            completed_jobs=$((completed_jobs + 1))
            log_ok "$name  (${elapsed}s)  [$completed_jobs/$total_jobs]"
            on_job_success "$name"
        else
            # ---- failure: retry or give up ----
            local retries
            retries=$(retries_get "$name")
            if [ "$retries" -lt "$MAX_RETRIES" ]; then
                retries=$((retries + 1))
                retries_set "$name" "$retries"
                state_set "$name" "PENDING"
                queue_push "$name"
                log_warn "$name failed (exit=$exit_code, ${elapsed}s) - retry $retries/$MAX_RETRIES queued"
            else
                state_set "$name" "FAILED"
                completed_jobs=$((completed_jobs + 1))
                failed_jobs=$((failed_jobs + 1))
                log_err "$name  (exit=$exit_code, ${elapsed}s, retries exhausted)"
            fi
        fi
    done
}

# ============================================================
#  Dependency chain: successful spinup to queue matching drainage
# ============================================================
on_job_success() {
    local name="$1"
    case "$name" in *_spinup) ;; *) return 0 ;; esac

    local core="${name%_spinup}"

    for f in "$CASE_DIR"/n[12]d*_drainage*.m; do
        [ -f "$f" ] || continue
        local dname
        dname=$(basename "${f%.m}")
        case_matches_workflow "$dname" || continue
        case "$dname" in "$core"_*_drainage*) ;;
            *) continue ;;
        esac
        if [ "$(state_get "$dname")" = "NONE" ]; then
            queue_push "$dname"
            retries_set "$dname" "0"
            state_set "$dname" "PENDING"
            total_jobs=$((total_jobs + 1))
            log "  Queued drainage: $dname  (total=$total_jobs)"
        fi
    done
}

# ============================================================
#  Summary
# ============================================================
print_summary() {
    echo ""
    echo "===================================================="
    echo "  JOB SUMMARY"
    echo "===================================================="
    local n_done=0 n_fail=0 n_kill=0 n_other=0
    for sf in "$STATE_DIR"/*.state; do
        [ -f "$sf" ] || continue
        local jname st
        jname=$(basename "${sf%.state}")
        st=$(cat "$sf")
        case "$st" in
            DONE)    n_done=$((n_done+1));   printf '  OK %s\n' "$jname" ;;
            FAILED)  n_fail=$((n_fail+1));   printf '  FAIL %s  -> %s/%s.log\n' "$jname" "$LOG_DIR" "$jname" ;;
            KILLED)  n_kill=$((n_kill+1));   printf '  KILLED %s\n' "$jname" ;;
            *)       n_other=$((n_other+1)); printf '  ? %s  (%s)\n' "$jname" "$st" ;;
        esac
    done
    echo "----------------------------------------------------"
    printf '  Done: %d   Failed: %d   Killed: %d\n' "$n_done" "$n_fail" "$n_kill"
    echo "===================================================="
}

# ============================================================
#  MAIN
# ============================================================
log "NEVIS Job Runner"
log "Config: parallel=$MAX_PARALLEL  retries=$MAX_RETRIES  timeout=${JOB_TIMEOUT}s  skip_spinup=$SKIP_SPINUP  case_dir=$CASE_DIR  workflow=$CASE_WORKFLOW"

# --- Prepare directories & clean previous state ---
mkdir -p "$LOG_DIR" "$STATE_DIR" "$DONE_DIR"
rm -f "$STATE_DIR"/*.state "$STATE_DIR"/*.pid "$STATE_DIR"/*.retries "$STATE_DIR"/*.start
: > "$QUEUE_FILE"

# --- Discover and enqueue scripts ---
if [ "$SKIP_SPINUP" = "true" ]; then
    for f in "$CASE_DIR"/n[12]d*.m; do
        [ -f "$f" ] || continue
        name=$(basename "${f%.m}")
        case_matches_workflow "$name" || continue
        case "$name" in *_spinup) continue ;; esac
        queue_push "$name"
        retries_set "$name" "0"
        state_set "$name" "PENDING"
        total_jobs=$((total_jobs + 1))
    done
    log "Skipping spinup - queued $total_jobs jobs directly"
else
    for f in "$CASE_DIR"/n[12]d*.m; do
        [ -f "$f" ] || continue
        name=$(basename "${f%.m}")
        case_matches_workflow "$name" || continue
        case "$name" in *_drainage*) continue ;; esac
        queue_push "$name"
        retries_set "$name" "0"
        state_set "$name" "PENDING"
        total_jobs=$((total_jobs + 1))
    done
    log "Queued $total_jobs initial jobs (drainage follows after spinup)"
fi

if [ "$total_jobs" -eq 0 ]; then
    log_err "No case configs found matching $CASE_DIR/n[12]d*.m - exiting."
    exit 1
fi
echo ""

# --- Main loop ---
while [ "$completed_jobs" -lt "$total_jobs" ]; do
    # Fill available slots from queue
    while [ "$(running_count)" -lt "$MAX_PARALLEL" ] && [ "$(queue_size)" -gt 0 ]; do
        if queue_pop; then
            launch_job "$REPLY"
        else
            break
        fi
    done

    sleep "$POLL_INTERVAL"
    reap_jobs

    # Stall detection
    if [ "$(running_count)" -eq 0 ] && [ "$(queue_size)" -eq 0 ] && [ "$completed_jobs" -lt "$total_jobs" ]; then
        log_warn "Stall: no running jobs, empty queue, $((total_jobs - completed_jobs)) jobs unfinished"
        log_warn "Possible cause: drainage waiting for a spinup that failed"
        break
    fi
done

print_summary

if [ "$failed_jobs" -gt 0 ]; then
    log_err "$failed_jobs jobs failed - check logs in $LOG_DIR/"
    exit 1
fi

log_ok "All $completed_jobs jobs completed successfully"
exit 0
