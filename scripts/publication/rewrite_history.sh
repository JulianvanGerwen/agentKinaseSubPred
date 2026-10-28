#!/usr/bin/env bash
# OPTIONAL rewrite of the git history before publication. READ docs/PUBLISHING.md section 5 FIRST.
#
# WARNING. This changes the hash of every rewritten commit and of every commit after it:
#   - FREEZE.json "harness_commit" (974d520, amendment A28) will no longer resolve;
#   - 406 results/runs/*/manifest.json "git_hash" values and the git_head in results/test/LEDGER.md will dangle;
#   - the tag targets in README.md (5fc7870, c3e2a14) must be looked up again in the commit map;
#   - anyone who already cloned the repository must clone it again.
# The script does not edit those files. It saves filter-repo's commit map (old hash -> new hash), which you should
# commit and publish so that the pre-registration can still be checked against the old hashes.
#
# Requirements:
#   - git-filter-repo, which is NOT installed in the project environment:  pip install git-filter-repo
#   - a FRESH CLONE made for this purpose. NEVER run this in the working project directory:
#         git clone --no-local <project-dir> ksp-rewrite && cd ksp-rewrite
#     filter-repo itself refuses to run in a repository that is not a fresh clone; do not pass --force to it.
#
# Usage (inside the fresh clone):
#   bash scripts/publication/rewrite_history.sh [--paths] [--manuscript] [--nfs] [--mailmap FILE] [--strip-session] --yes
#     --paths          remove the 13 condition tables of the provided compendium from ALL commits
#                      (11 paths at HEAD + f6b_conditions_lam10/lam100.tsv in history; docs/PUBLISHING.md section 3)
#     --manuscript     remove manuscript/ (the paper draft, published separately) from ALL commits
#     --nfs            replace the internal file-server address and share path with '<internal-file-server>:<share>'
#                      in all blobs (one review file, docs/reviews/r1_ai_failure_modes.md, in 103 commits)
#     --mailmap FILE   rewrite author/committer identities with a git mailmap file that you write yourself, e.g.
#                      'Public Name <public address> <address used in the commits>'
#     --strip-session  remove the 'Claude-Session: <URL>' trailer lines from commit messages
#     --yes            required; confirms that you read the warnings above
set -euo pipefail

DO_PATHS=0; DO_MS=0; DO_NFS=0; MAILMAP=""; DO_SESSION=0; YES=0
while [ $# -gt 0 ]; do
  case "$1" in
    --paths) DO_PATHS=1 ;;
    --manuscript) DO_MS=1 ;;
    --nfs) DO_NFS=1 ;;
    --mailmap) MAILMAP="${2:?--mailmap needs a file}"; shift ;;
    --strip-session) DO_SESSION=1 ;;
    --yes) YES=1 ;;
    -h|--help) sed -n '2,33p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
  shift
done
[ "$YES" = 1 ] || { echo "refusing: pass --yes after reading the warnings (bash $0 --help)" >&2; exit 2; }
[ $((DO_PATHS + DO_MS + DO_NFS + DO_SESSION)) -gt 0 ] || [ -n "$MAILMAP" ] || { echo "nothing to do" >&2; exit 2; }
command -v git-filter-repo >/dev/null 2>&1 || git filter-repo --version >/dev/null 2>&1 || {
  echo "git-filter-repo is not installed: pip install git-filter-repo" >&2; exit 1; }

TOP="$(git rev-parse --show-toplevel)"
cd "$TOP"
# A working project directory has the restricted inputs next to the code; a fresh clone does not.
# (sealed/ is not checked: a clone contains the tracked sealed/kinase_label_counts.tsv.)
for d in data external cache; do
  if [ -e "$d" ]; then
    echo "refusing: '$d/' exists here, so this looks like the working project, not a fresh clone" >&2; exit 1
  fi
done
if [ -n "$(git status --porcelain)" ]; then echo "refusing: the working tree is not clean" >&2; exit 1; fi

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
ARGS=()

if [ "$DO_PATHS" = 1 ] || [ "$DO_MS" = 1 ]; then
  : > "$WORK/paths.txt"
  [ "$DO_MS" = 1 ] && echo "manuscript/" >> "$WORK/paths.txt"
fi
if [ "$DO_PATHS" = 1 ]; then
  cat >> "$WORK/paths.txt" <<'EOF'
results/data_profile/compendium_per_condition_coverage.tsv
results/data_profile/kinhib_hijazi_conditions.tsv
results/data_profile/kinhib_condition_matches.tsv
results/data_profile/compendium_conditions_by_pmid.tsv
results/features/phospho/conditions_meta.tsv
results/features/phospho/f6b_conditions.tsv
results/features/phospho/f6b_conditions_clean_bench.tsv
results/features/phospho/f6b_conditions_f6_clean.tsv
results/features/phospho/f6b_conditions_lam10.tsv
results/features/phospho/f6b_conditions_lam100.tsv
results/features/phospho/f6b_unmapped_or_dropped_drug_conditions.tsv
results/features/phospho/label_source_studies.tsv
results/features/phospho/benchmarkin_studies.tsv
EOF
fi
if [ "$DO_PATHS" = 1 ] || [ "$DO_MS" = 1 ]; then
  ARGS+=(--invert-paths --paths-from-file "$WORK/paths.txt")
fi

if [ "$DO_NFS" = 1 ]; then
  # An RFC 1918 10.x.y.z address directly followed by ':/<share path>'. The pattern is generic on purpose, so that
  # the address itself is not written into the repository.
  printf '%s\n' 'regex:(?<![0-9.])10\.[0-9]{1,3}\.[0-9]{1,3}\.[0-9]{1,3}:/[^\s`'"'"'")]*==><internal-file-server>:<share>' \
    > "$WORK/replace.txt"
  ARGS+=(--replace-text "$WORK/replace.txt")
fi

if [ "$DO_SESSION" = 1 ]; then
  printf '%s\n' 'regex:\n?Claude-Session: https://claude\.ai/code/session_[A-Za-z0-9]+==>' > "$WORK/msg.txt"
  ARGS+=(--replace-message "$WORK/msg.txt")
fi

if [ -n "$MAILMAP" ]; then
  ARGS+=(--mailmap "$MAILMAP")
fi

echo "git filter-repo ${ARGS[*]}"
git filter-repo "${ARGS[@]}"

mkdir -p docs
cp .git/filter-repo/commit-map docs/history_commit_map.tsv
echo
echo "Done. The old -> new commit map is in docs/history_commit_map.tsv (not committed yet)."
echo "Next steps (docs/PUBLISHING.md section 5):"
echo "  1. look up the new hashes of 974d520, 5fc7870 and c3e2a14 in docs/history_commit_map.tsv;"
echo "  2. git add docs/history_commit_map.tsv && git commit -m 'Add commit map of the history rewrite'"
echo "  3. create the tags on the NEW hashes, then push to a NEW remote (git remote add origin <url>)."
