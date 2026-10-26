#!/bin/bash
# Sequential queue of 5-fold configs. Lines: name|threads|stage|sets   (waits for a pid first if WAIT_PID set)
cd /workspace
[ -n "$WAIT_PID" ] && while kill -0 $WAIT_PID 2>/dev/null; do sleep 20; done
while IFS='|' read -r name thr stage sets; do
  [ -z "$name" ] && continue; [[ "$name" == \#* ]] && continue
  echo "=== $(date +%H:%M:%S) start $name"
  scripts/models/r2_cv5.sh "$name" "$thr" "$stage" "$sets"
done < "$1"
