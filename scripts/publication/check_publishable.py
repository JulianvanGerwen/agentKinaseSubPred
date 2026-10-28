#!/usr/bin/env python3
"""Pre-publication hygiene scan of the files that a push would publish.

Usage (from the repository root):
    python -I scripts/publication/check_publishable.py                 # walk the tree, apply .gitignore
    python -I scripts/publication/check_publishable.py --git           # use `git ls-files -co --exclude-standard`
    python -I scripts/publication/check_publishable.py --files-from L  # scan the paths listed in file L
    python -I scripts/publication/check_publishable.py --check-ignore P [P ...]   # which paths .gitignore matches
    python -I scripts/publication/check_publishable.py --check-ignore-from L      # same, paths listed in file L

What it looks for (ERROR = exit code 1):
  - e-mail addresses; the unix user and group names of whoever runs it; the machine's host name
  - absolute home paths (/home/..., /Users/..., /root/...) and agent scratchpad paths (/tmp/claude-...)
  - private IPv4 addresses (10/8, 172.16/12, 192.168/16); session URLs (claude.ai/code/session_...)
  - secret patterns (cloud and API keys, private-key blocks, password=/token= assignments)
  - files larger than 50 MB; file names of restricted inputs (PhosphoSitePlus tables, compendium CSVs, arrays)
WARN (reported, exit code unaffected): absolute /workspace paths; files larger than 5 MB; inside
results/runs/*/manifest.json (run records, kept unedited on purpose) the host field, the identity names, home
paths and scratchpad paths.

No third-party dependency. The .gitignore matcher implements the usual git rules (anchoring, '**', negation,
directory-only patterns, an ignored directory hides everything below it); `pathspec` is used instead if installed.
"""
from __future__ import annotations

import argparse
import getpass
import grp
import os
import re
import socket
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BIG_ERROR = 50 * 1024 * 1024
BIG_WARN = 5 * 1024 * 1024
SKIP_DIRS = {".git", ".venv", "__pycache__", ".pytest_cache"}
TEXT_LIMIT = 20 * 1024 * 1024           # larger files are checked for size only
SELF = Path(__file__).resolve()

# ------------------------------------------------------------------------------------------- .gitignore matching


def _glob_to_re(pat: str) -> str:
    i, out = 0, []
    while i < len(pat):
        c = pat[i]
        if c == "*":
            if pat[i:i + 2] == "**":
                # '**/' = zero or more directories; trailing '/**' = everything inside
                if pat[i:i + 3] == "**/":
                    out.append(r"(?:.*/)?")
                    i += 3
                    continue
                out.append(r".*")
                i += 2
                continue
            out.append(r"[^/]*")
        elif c == "?":
            out.append(r"[^/]")
        elif c == "[":
            j = pat.find("]", i + 1)
            if j == -1:
                out.append(r"\[")
            else:
                body = pat[i + 1:j].replace("\\", "\\\\")
                if body.startswith("!"):
                    body = "^" + body[1:]
                out.append(f"[{body}]")
                i = j
        elif c == "\\" and i + 1 < len(pat):
            out.append(re.escape(pat[i + 1]))
            i += 1
        else:
            out.append(re.escape(c))
        i += 1
    return "".join(out)


class GitIgnore:
    def __init__(self, text: str):
        self.rules = []
        for raw in text.splitlines():
            line = raw.rstrip()
            if not line or line.startswith("#"):
                continue
            neg = line.startswith("!")
            if neg:
                line = line[1:]
            dir_only = line.endswith("/")
            line = line.rstrip("/")
            anchored = "/" in line
            line = line.lstrip("/")
            body = _glob_to_re(line)
            rx = re.compile(("^" if anchored else r"^(?:.*/)?") + body + "$")
            self.rules.append((rx, neg, dir_only, raw))

    def _match(self, rel: str, is_dir: bool):
        hit = None
        for rx, neg, dir_only, raw in self.rules:
            if dir_only and not is_dir:
                continue
            if rx.match(rel):
                hit = None if neg else raw
        return hit

    def ignored(self, rel: str, is_dir: bool = False):
        """Return the matching rule (str) if rel is ignored, else None."""
        parts = rel.split("/")
        for k in range(1, len(parts)):
            r = self._match("/".join(parts[:k]), True)
            if r:
                return r
        return self._match(rel, is_dir)


def load_ignore(root: Path):
    f = root / ".gitignore"
    text = f.read_text() if f.exists() else ""
    try:
        import pathspec  # type: ignore

        spec = pathspec.GitIgnoreSpec.from_lines(text.splitlines())

        class _PS:
            def ignored(self, rel, is_dir=False):
                return "pathspec" if spec.match_file(rel + ("/" if is_dir else "")) else None
        return _PS()
    except ImportError:
        return GitIgnore(text)


def walk_files(root: Path, ign) -> list[str]:
    out = []
    for dp, dns, fns in os.walk(root):
        rd = Path(dp).relative_to(root).as_posix()
        rd = "" if rd == "." else rd
        keep = []
        for d in dns:
            rel = f"{rd}/{d}" if rd else d
            if d in SKIP_DIRS or ign.ignored(rel, True):
                continue
            keep.append(d)
        dns[:] = keep
        for fn in fns:
            rel = f"{rd}/{fn}" if rd else fn
            if not ign.ignored(rel, False):
                out.append(rel)
    return sorted(out)


# ------------------------------------------------------------------------------------------------- patterns


def identity_terms() -> list[str]:
    terms = set()
    try:
        terms.add(getpass.getuser())
    except Exception:
        pass
    try:
        terms.add(grp.getgrgid(os.getgid()).gr_name)
    except Exception:
        pass
    h = socket.gethostname()
    terms.update({h, h.split(".")[0]})
    return sorted(t for t in terms if t and len(t) >= 4 and t not in {"root", "user", "localhost"})


EMAIL = re.compile(r"(?<![\w.+-])[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)*\.[A-Za-z]{2,}\b")
EMAIL_OK = re.compile(r"^(?:noreply@anthropic\.com|git@github\.com|[^@]+@users\.noreply\.github\.com)$")
HOME = re.compile(r"(?<![\w.~-])/(?:home|Users|root)/(?!\.\.\.)[A-Za-z0-9._-]+")
SCRATCH = re.compile(r"/tmp/claude-[0-9]+[^\s`'\")]*")
WORKSPACE = re.compile(r"(?<![\w.~-])/workspace\b")
PRIV_IP = re.compile(r"(?<![\d.])(?:10\.\d{1,3}|172\.(?:1[6-9]|2\d|3[01])|192\.168)\.\d{1,3}\.\d{1,3}(?![\d.])")
SESSION = re.compile(r"claude\.ai/code/session_[A-Za-z0-9]+")
SECRETS = [
    ("aws-key", re.compile(r"\bAKIA[0-9A-Z]{16}\b")),
    ("github-token", re.compile(r"\b(?:ghp|gho|ghu|ghs|ghr)_[A-Za-z0-9]{30,}\b|\bgithub_pat_[A-Za-z0-9_]{20,}")),
    ("anthropic/openai-key", re.compile(r"\bsk-(?:ant-)?[A-Za-z0-9_-]{24,}")),
    ("huggingface-token", re.compile(r"\bhf_[A-Za-z0-9]{30,}\b")),
    ("slack-token", re.compile(r"\bxox[abprs]-[A-Za-z0-9-]{10,}")),
    ("private-key", re.compile(r"-----BEGIN (?:[A-Z]+ )?PRIVATE KEY-----")),
    ("assignment", re.compile(r"(?i)\b(?:password|passwd|api[_-]?key|secret|access[_-]?token|auth[_-]?token)\s*[=:]\s*['\"][^'\"\s]{8,}['\"]")),
    ("bearer", re.compile(r"\bBearer\s+[A-Za-z0-9._~+/-]{20,}=*")),
]
RESTRICTED_NAMES = [
    re.compile(r"(?i)(?:^|/)(?:Kinase_Substrate_Dataset|Phosphorylation_site_dataset|Regulatory_sites|"
               r"Disease-associated_sites|Phosphosite_PTM_seq|Phosphosite_seq)(?:\.gz|\.txt|\.tsv)?$"),
    re.compile(r"(?i)(?:^|/)0[1-3]_phospho\w*\.csv$"),
    re.compile(r"(?i)\.(?:parquet|npy|npz|pt|pth|ckpt|safetensors|h5|pkl|feather|rda|Rds)$"),
    re.compile(r"(?:^|/)sealed/(?!kinase_label_counts\.tsv$)"),
    re.compile(r"(?:^|/)(?:data|external|cache)/"),
]
MANIFEST = re.compile(r"^results/runs/[^/]+/manifest\.json$")
HOST_FIELD = re.compile(r'"host"\s*:\s*"[^"]+"')


def is_text(p: Path) -> bool:
    try:
        with open(p, "rb") as fh:
            return b"\0" not in fh.read(8192)
    except OSError:
        return False


def scan(root: Path, files: list[str]):
    errors, warns = [], []
    ids = identity_terms()
    id_rx = re.compile(r"(?i)(?<![A-Za-z0-9])(?:" + "|".join(map(re.escape, ids)) + r")(?![A-Za-z0-9])") if ids else None
    for rel in files:
        p = root / rel
        if not p.is_file() or p.resolve() == SELF:
            continue
        size = p.stat().st_size
        if size > BIG_ERROR:
            errors.append((rel, 0, "file>50MB", f"{size / 1e6:.1f} MB"))
        elif size > BIG_WARN:
            warns.append((rel, 0, "file>5MB", f"{size / 1e6:.1f} MB"))
        for rx in RESTRICTED_NAMES:
            if rx.search(rel):
                errors.append((rel, 0, "restricted-filename", rx.pattern[:60]))
                break
        if size > TEXT_LIMIT or not is_text(p):
            continue
        text = p.read_text(errors="replace")
        manifest = bool(MANIFEST.match(rel))
        n_ws = 0
        for ln, line in enumerate(text.splitlines(), 1):
            for m in EMAIL.finditer(line):
                if not EMAIL_OK.match(m.group(0)):
                    errors.append((rel, ln, "email", m.group(0)))
            if id_rx:
                for m in id_rx.finditer(line):
                    (warns if manifest else errors).append((rel, ln, "user/group/host name", m.group(0)))
            for m in HOME.finditer(line):
                (warns if manifest else errors).append((rel, ln, "home-path", m.group(0)))
            for m in SCRATCH.finditer(line):
                (warns if manifest else errors).append((rel, ln, "scratchpad-path", m.group(0)[:60]))
            for m in PRIV_IP.finditer(line):
                errors.append((rel, ln, "private-ip", m.group(0)))
            for m in SESSION.finditer(line):
                errors.append((rel, ln, "session-url", m.group(0)[:40]))
            for name, rx in SECRETS:
                for m in rx.finditer(line):
                    errors.append((rel, ln, f"secret:{name}", m.group(0)[:12] + "..."))
            n_ws += len(WORKSPACE.findall(line))
            if manifest and HOST_FIELD.search(line):
                warns.append((rel, ln, "host field in run manifest", "(value not shown)"))
        if n_ws:
            warns.append((rel, 0, "absolute /workspace path", f"{n_ws} occurrence(s)"))
    return errors, warns


def git_files(root: Path) -> list[str]:
    r = subprocess.run(["git", "ls-files", "-co", "--exclude-standard", "-z"], cwd=root, capture_output=True)
    if r.returncode != 0:
        sys.exit("git ls-files failed: " + r.stderr.decode(errors="replace"))
    return sorted(x for x in r.stdout.decode().split("\0") if x)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=str(REPO))
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--git", action="store_true", help="list files with git (tracked + untracked, not ignored)")
    src.add_argument("--files-from", help="file with one repository-relative path per line")
    src.add_argument("--check-ignore", nargs="+", metavar="PATH", help="print the .gitignore rule matching each path")
    src.add_argument("--check-ignore-from", metavar="FILE", help="like --check-ignore, paths read from FILE")
    ap.add_argument("--max-show", type=int, default=200, help="findings printed per severity")
    a = ap.parse_args(argv)
    root = Path(a.root).resolve()
    ign = load_ignore(root)

    if a.check_ignore or a.check_ignore_from:
        paths = a.check_ignore or [x.strip() for x in Path(a.check_ignore_from).read_text().splitlines() if x.strip()]
        n = 0
        for rel in paths:
            r = ign.ignored(rel.rstrip("/"), rel.endswith("/") or (root / rel).is_dir())
            if r:
                n += 1
                print(f"IGNORED\t{rel}\t{r}")
        print(f"# {n} of {len(paths)} paths match .gitignore", file=sys.stderr)
        return 0

    if a.git:
        files = git_files(root)
    elif a.files_from:
        files = [x.strip() for x in Path(a.files_from).read_text().splitlines() if x.strip()]
    else:
        files = walk_files(root, ign)
    errors, warns = scan(root, files)
    total = sum((root / f).stat().st_size for f in files if (root / f).is_file())
    print(f"scanned {len(files)} files, {total / 1e6:.1f} MB")
    for label, rows in (("ERROR", errors), ("WARN", warns)):
        print(f"\n{label}: {len(rows)}")
        for rel, ln, kind, val in rows[:a.max_show]:
            print(f"  {label}\t{kind}\t{rel}:{ln}\t{val}")
        if len(rows) > a.max_show:
            print(f"  ... {len(rows) - a.max_show} more")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
