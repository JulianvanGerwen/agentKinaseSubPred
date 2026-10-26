"""Tiny job queue for the regime runs: launches shell commands from a queue file (one per line, '#' = comment) while
(#running < limit in <queue>.limit) and (total RSS of this user's regime_* processes < --max-rss GB).
Usage: regime_queue.py --queue logs/regime/xfit.q --max-rss 55"""
import argparse, os, subprocess, time
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--queue", required=True)
ap.add_argument("--max-rss", type=float, default=55.0)
a = ap.parse_args()
q = Path(a.queue)
lim_p = q.with_suffix(".limit")
done_p = q.with_suffix(".done")
running = {}


def rss_gb():
    out = subprocess.run("ps -u $(id -u) -o rss=,args= | grep regime_ | grep -v grep | awk '{s+=$1} END {print s}'",
                         shell=True, capture_output=True, text=True).stdout.strip()
    return float(out or 0) / 1e6


lines = [l.strip() for l in q.read_text().splitlines() if l.strip() and not l.startswith("#")]
i = 0
while i < len(lines) or running:
    for p in list(running):
        if p.poll() is not None:
            with open(done_p, "a") as fh:
                fh.write(f"{time.strftime('%H:%M:%S')} rc={p.returncode} {running[p]}\n")
            del running[p]
    limit = int(lim_p.read_text().strip()) if lim_p.exists() else 3
    if i < len(lines) and len(running) < limit and rss_gb() < a.max_rss:
        cmd = lines[i]
        i += 1
        running[subprocess.Popen(cmd, shell=True, cwd="/workspace")] = cmd
        with open(done_p, "a") as fh:
            fh.write(f"{time.strftime('%H:%M:%S')} START {cmd}\n")
        time.sleep(20)
        continue
    time.sleep(15)
