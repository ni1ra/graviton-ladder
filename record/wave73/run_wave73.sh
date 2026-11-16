#!/usr/bin/env bash
# Reruns the quick, deterministic checks of record/wave73 and fails if any of them changes its verdict.
#
#   record/wave73/run_wave73.sh                 every check in wave73.manifest
#   ONLY='even-(dict|const)' ...run_wave73.sh   only the checks whose id matches the regex
#   JOBS=4 ...run_wave73.sh                     slots (default: the number of cores, at most 8). A check holds its `threads`
#                                               slots (default 1) and gets that many BLAS threads. JOBS=1: one check at a time,
#                                               each on every core
#   SECTION_BUDGET=600 ...run_wave73.sh         fail, naming what still ran, if the section takes longer than this many seconds
#   SLOW_CHECKS=1 ...run_wave73.sh              also the checks marked `slow` in the manifest (kept out of the default run for time)
#   KEEP=1 ...run_wave73.sh                     keep the temporary directory and print where it is
#   MANIFEST=/path ...run_wave73.sh             another manifest (used to prove that a wrong one turns red)
#   PYTHON=python3                              the interpreter (python-flint is needed by the flint checks)
#
# What it does, per check of the manifest (format documented at the top of wave73.manifest):
#   1. copies record/wave73 into a temporary directory, so nothing is ever written into the repository
#      (several historical scripts write out/ and logs/ next to themselves);
#   2. runs the command from the folder the manifest names, output to a file in the temporary directory;
#   3. compares the exit code with the manifest. Several scripts exit 1 BY DESIGN because one
#      pre-registered line failed and the report says so; that is part of the manifest, not a bug;
#   4. requires every `need` line (a regular expression, matched against whole output lines);
#   5. requires each `same` / `contains` stored output to be reproduced (analyses of stored raw rows)
#   6. requires every `count` to match exactly (so a check that silently disappears turns CI red);
#   7. treats every output line that mentions FAIL, FAILED, MISSED, MISMATCH, Traceback or Error as a
#      failure unless it is matched by one of the check's own `need` lines. So a script that is supposed to
#      report one known failed line and starts reporting two, or zero, is red.
# It also fails if a folder of record/wave73 is neither checked nor named by a `skip` line, so a new folder
# cannot be added without a decision about whether it is rerun.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANIFEST="${MANIFEST:-$HERE/wave73.manifest}"
PY="${PYTHON:-python3}"
if [ -z "${JOBS:-}" ]; then
  JOBS="$(nproc 2>/dev/null || echo 1)"; [ "$JOBS" -gt 8 ] && JOBS=8
  # a container limited with `docker run --cpus N` still reports every CPU of the host to nproc
  quota="$(awk '$1 != "max" && $2 > 0 { print int(($1 + $2 - 1) / $2) }' /sys/fs/cgroup/cpu.max 2>/dev/null || true)"
  if [ -n "$quota" ] && [ "$quota" -ge 1 ] && [ "$quota" -lt "$JOBS" ]; then JOBS="$quota"; fi
fi
ONLY="${ONLY:-}"
CHECK_TIMEOUT="${CHECK_TIMEOUT:-1800}"
BAD='(^|[^A-Za-z])(FAIL|FAILED|MISSED|MISMATCH)([^A-Za-z]|$)|Traceback|Error'

NL=$'\n'

[ -f "$MANIFEST" ] || { echo "run_wave73: manifest not found: $MANIFEST" >&2; exit 2; }

TMP="$(mktemp -d "${TMPDIR:-/tmp}/wave73.XXXXXX")"
if [ -n "${KEEP:-}" ]; then echo "run_wave73: temporary directory kept: $TMP"; else trap 'rm -rf "$TMP"' EXIT; fi
cp -R "$HERE" "$TMP/wave73"
mkdir -p "$TMP/res" "$TMP/art"

export PYTHONDONTWRITEBYTECODE=1 PYTHONUTF8=1 PYTHONIOENCODING=utf-8 LC_ALL="${LC_ALL:-C.UTF-8}"

# ---- parse the manifest ------------------------------------------------------------------------------
n=0
declare -a ID AFTER THREADS SLOW PID CWD CMD RC NEED COUNT SAME STARTS
declare -a SKIPDIRS=()
cur=-1
lineno=0
while IFS= read -r line || [ -n "$line" ]; do
  lineno=$((lineno + 1))
  line="${line%$'\r'}"
  case "$line" in ''|'#'*) continue ;; esac
  kw="${line%% *}"; rest="${line#* }"
  [ "$kw" = "$line" ] && rest=""
  case "$kw" in
    check) cur=$n; n=$((n + 1)); ID[$cur]="$rest"; AFTER[$cur]=""; THREADS[$cur]=1; SLOW[$cur]=0; RC[$cur]=0; NEED[$cur]=""; COUNT[$cur]=""; SAME[$cur]=""; STARTS[$cur]=""; CWD[$cur]=""; CMD[$cur]="" ;;
    after) AFTER[$cur]="$rest" ;;
    threads) THREADS[$cur]="$rest" ;;
    slow)  SLOW[$cur]=1 ;;
    run)   CWD[$cur]="${rest%% *}"; CMD[$cur]="${rest#* }" ;;
    rc)    RC[$cur]="$rest" ;;
    need)  NEED[$cur]+="$rest"$'\n' ;;
    count) COUNT[$cur]+="$rest"$'\n' ;;
    same)  SAME[$cur]+="$rest"$'\n' ;;
    contains) STARTS[$cur]+="$rest"$'\n' ;;
    skip)  SKIPDIRS+=("${rest%% *}") ;;
    *) echo "run_wave73: $MANIFEST:$lineno: unknown keyword '$kw'" >&2; exit 2 ;;
  esac
done < "$MANIFEST"

# ---- every folder must be checked or explicitly skipped --------------------------------------------
uncovered=0
for d in "$HERE"/*/; do
  d="$(basename "$d")"
  covered=0
  for ((i = 0; i < n; i++)); do [ "${CWD[$i]%%/*}" = "$d" ] && covered=1; done
  for s in "${SKIPDIRS[@]:-}"; do [ "${s%%/*}" = "$d" ] && covered=1; done
  if [ $covered = 0 ]; then echo "FAIL coverage: folder $d is neither checked nor skipped in $(basename "$MANIFEST")"; uncovered=1; fi
done

# ---- run -------------------------------------------------------------------------------------------
run_one() {
  local i=$1 id="${ID[$1]}" art="$TMP/art/${ID[$1]}"
  mkdir -p "$art" "$TMP/wave73/${CWD[$i]%%/*}/out" "$TMP/wave73/${CWD[$i]%%/*}/logs"
  local cmd="${CMD[$i]//\{OUT\}/$art}"
  local t0=$SECONDS
  local th=""; [ "$JOBS" -gt 1 ] && th="${THREADS[$i]}"
  ( cd "$TMP/wave73/${CWD[$i]}" && { [ -z "$th" ] || export OMP_NUM_THREADS="$th" OPENBLAS_NUM_THREADS="$th" MKL_NUM_THREADS="$th"; } && timeout "$CHECK_TIMEOUT" "$PY" -B $cmd ) > "$TMP/res/$id.out" 2>&1
  echo $? > "$TMP/res/$id.rc"
  echo $((SECONDS - t0)) > "$TMP/res/$id.sec"
}

failed=0; ran=0; sum_checks=0; sum_plants=0; t_start=$SECONDS

evaluate() {
  local i=$1 id="${ID[$1]}" out="$TMP/res/${ID[$1]}.out"
  local rc sec why="" lines
  rc="$(cat "$TMP/res/$id.rc")"; sec="$(cat "$TMP/res/$id.sec")"
  [ "$rc" = "${RC[$i]}" ] || why+=$'\n    exit code '"$rc"', manifest says '"${RC[$i]}"
  while IFS= read -r re; do
    [ -z "$re" ] && continue
    grep -Eq -- "$re" "$out" || why+=$'\n    required line missing: '"$re"
  done <<< "${NEED[$i]}"
  local nchk=0 npl=0
  while IFS= read -r spec; do
    [ -z "$spec" ] && continue
    # count LABEL N REGEX [ ::: EXCLUDE-REGEX ]
    local label="${spec%% *}"; spec="${spec#* }"; local want="${spec%% *}"; spec="${spec#* }"
    local re="${spec%% ::: *}" ex=""; [ "$re" != "$spec" ] && ex="${spec#* ::: }"
    local got
    if [ -n "$ex" ]; then got="$(grep -E -- "$re" "$out" | grep -Evc -- "$ex")"; else got="$(grep -Ec -- "$re" "$out")"; fi
    [ "$got" = "$want" ] || why+=$'\n    count '"$label"': '"$got"' lines match /'"$re"'/, manifest says '"$want"
    case "$label" in checks) nchk=$((nchk + want)) ;; plants) npl=$((npl + want)) ;; esac
  done <<< "${COUNT[$i]}"
  # stored outputs that the run must reproduce: `same` = identical, `contains` = every line of the stored file is
  # reproduced in order (an addendum the script gained after the stored run is allowed). Carriage returns of Windows-written files are ignored.
  local exempt="$TMP/res/$id.exempt" norm="$TMP/res/$id.norm" ref="$TMP/res/$id.ref"; : > "$exempt"
  tr -d '\015' < "$out" > "$norm"
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    if [ ! -f "$HERE/$f" ]; then why+="${NL}    stored file missing: $f"; continue; fi
    tr -d '\015' < "$HERE/$f" > "$ref"; cat "$ref" >> "$exempt"
    cmp -s "$norm" "$ref" || why+="${NL}    output differs from stored $f${NL}$(diff "$ref" "$norm" | head -n 6 | cut -c1-160 | sed 's/^/      /')"
  done <<< "${SAME[$i]}"
  while IFS= read -r f; do
    [ -z "$f" ] && continue
    if [ ! -f "$HERE/$f" ]; then why+="${NL}    stored file missing: $f"; continue; fi
    tr -d '\015' < "$HERE/$f" > "$ref"; cat "$ref" >> "$exempt"
    # no line of the stored file may be missing or changed; lines the current script adds are allowed
    [ "$(diff "$ref" "$norm" | grep -c '^<')" = 0 ] || why+="${NL}    output lacks or changes lines of stored $f${NL}$(diff "$ref" "$norm" | grep '^<' | head -n 6 | cut -c1-160 | sed 's/^/      /')"
  done <<< "${STARTS[$i]}"
  # every line that looks like a failure must be one the manifest expects
  while IFS= read -r bad; do
    [ -z "$bad" ] && continue
    local ok=0
    grep -Fxq -- "$bad" "$exempt" && ok=1
    while IFS= read -r re; do [ -n "$re" ] && printf '%s\n' "$bad" | grep -Eq -- "$re" && { ok=1; break; }; done <<< "${NEED[$i]}"
    [ $ok = 1 ] || why+=$'\n    unexpected failure-like line: '"${bad:0:200}"
  done < <(grep -E -- "$BAD" "$out" || true)
  ran=$((ran + 1)); sum_checks=$((sum_checks + nchk)); sum_plants=$((sum_plants + npl))
  if [ -z "$why" ]; then
    printf 'ok    %-22s %4ss  checks=%s plants=%s  [section %ss]\n' "$id" "$sec" "$nchk" "$npl" "$((SECONDS - t_start))"
  else
    failed=$((failed + 1))
    printf 'FAIL  %-22s %4ss  [section %ss]%s\n' "$id" "$sec" "$((SECONDS - t_start))" "$why"
    echo "    --- last lines of its output ---"; tail -n 12 "$out" | cut -c1-200 | sed 's/^/    | /'
  fi
}

# Scheduler: up to JOBS checks at once, in manifest order (longest first keeps the machine busy), but a check
# with an `after` line starts only once the checks it names have finished.
pending=()
skipped_slow=0
wanted() {   # is check $1 part of this run?
  [ -n "$ONLY" ] && ! [[ "${ID[$1]}" =~ $ONLY ]] && return 1
  if [ "${SLOW[$1]}" = 1 ] && [ -z "${SLOW_CHECKS:-}" ]; then skipped_slow=$((skipped_slow + 1)); return 1; fi
  return 0
}
for ((i = 0; i < n; i++)); do wanted "$i" && pending+=("$i"); done
selected=" ${pending[*]:-} "
ready() {   # all selected dependencies of check $1 have finished
  local dep
  for dep in ${AFTER[$1]//,/ }; do
    for ((j = 0; j < n; j++)); do
      if [ "${ID[$j]}" = "$dep" ] && [[ "$selected" == *" $j "* ]] && [ ! -f "$TMP/res/$dep.sec" ]; then return 1; fi
    done
  done
  return 0
}
active=()   # checks that have been started and not yet seen finished
used=0      # slots they hold (a check holds `threads` of the JOBS slots)
refresh() {   # a check that has finished is judged and printed at once, so a timeout shows where the time went
  local keep=() a
  used=0
  for a in ${active[@]+"${active[@]}"}; do
    if [ -f "$TMP/res/${ID[$a]}.sec" ]; then evaluate "$a"; continue; fi
    keep+=("$a"); used=$((used + THREADS[a]))
  done
  active=(${keep[@]+"${keep[@]}"})
}
over_budget() {   # SECTION_BUDGET seconds: fail with names instead of being cancelled silently by the caller
  [ -n "${SECTION_BUDGET:-}" ] && [ $((SECONDS - t_start)) -ge "$SECTION_BUDGET" ] || return 1
  echo "FAIL budget: the section passed its budget of ${SECTION_BUDGET} s ($((SECONDS - t_start)) s in)"
  local a names=""
  for a in ${active[@]+"${active[@]}"}; do names+=" ${ID[$a]}"; pkill -P "${PID[$a]}" 2>/dev/null; kill "${PID[$a]}" 2>/dev/null; done
  echo "  still running:${names:- none}"
  names=""; for a in ${pending[@]+"${pending[@]}"}; do names+=" ${ID[$a]}"; done
  echo "  not started:${names:- none}"
  return 0
}
while :; do
  refresh
  [ ${#pending[@]} -eq 0 ] && [ ${#active[@]} -eq 0 ] && break
  if over_budget; then exit 1; fi
  launched=0
  for k in "${!pending[@]}"; do
    i="${pending[$k]}"
    # a check that needs more slots than JOBS still runs, alone
    if [ $((used + THREADS[i])) -le "$JOBS" ] || [ "$used" -eq 0 ]; then
      if ready "$i"; then
        run_one "$i" &
        PID[$i]=$!
        active+=("$i"); launched=1
        unset 'pending[k]'; pending=(${pending[@]+"${pending[@]}"})
        break
      fi
    fi
  done
  if [ $launched = 0 ]; then
    [ ${#active[@]} -gt 0 ] || { echo "run_wave73: circular or unknown 'after' in the manifest" >&2; exit 2; }
    sleep 0.5   # polling, not wait -n, so the budget is checked even while one check hangs
  fi
done
[ "$skipped_slow" -gt 0 ] && echo "wave73: $skipped_slow slow checks listed in the manifest were not run (SLOW_CHECKS=1 runs them)"
echo "wave73: $ran checks run in $((SECONDS - t_start)) s; $sum_checks check lines and $sum_plants planted errors as the manifest requires; $failed failed"
if [ "$failed" -ne 0 ] || [ "$uncovered" -ne 0 ]; then exit 1; fi
if [ "$ran" -eq 0 ]; then echo "wave73: no check ran" >&2; exit 1; fi
exit 0
