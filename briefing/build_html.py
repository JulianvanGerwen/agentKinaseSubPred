"""Build a self-contained briefing.html from briefing.md.

The `markdown` package is not installed in the project environment, so this
script contains a small converter for the subset of Markdown that
briefing.md uses: headings, paragraphs, bold, italic, inline code,
ordered and unordered lists, block quotes (nested Markdown inside),
pipe tables, horizontal rules and images. Images are embedded as base64,
and an image followed by a paragraph that starts with "**Figure" becomes a
<figure> with a caption. No external scripts, fonts or stylesheets.

Run from anywhere:  python -I briefing/build_html.py
"""
import base64
import html
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "briefing.md")
OUT = os.path.join(HERE, "briefing.html")

CSS = """
:root {
  color-scheme: light;
  --bg: #fcfcfb; --surface: #ffffff; --ink: #1a1a19; --ink2: #52514e;
  --rule: #dddcd7; --accent: #2a78d6; --flag-bg: #fdf0ea; --flag-edge: #eb6834;
  --box-bg: #f4f3f0; --code-bg: #efeeea;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --bg: #1a1a19; --surface: #242422; --ink: #f0efec; --ink2: #c3c2b7;
    --rule: #3a3a37; --accent: #6da7ec; --flag-bg: #3a2a22; --flag-edge: #d95926;
    --box-bg: #262624; --code-bg: #2e2e2b;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --bg: #1a1a19; --surface: #242422; --ink: #f0efec; --ink2: #c3c2b7;
  --rule: #3a3a37; --accent: #6da7ec; --flag-bg: #3a2a22; --flag-edge: #d95926;
  --box-bg: #262624; --code-bg: #2e2e2b;
}
* { box-sizing: border-box; }
html { -webkit-text-size-adjust: 100%; }
body {
  margin: 0; padding: 24px 16px 64px; background: var(--bg); color: var(--ink);
  font: 17px/1.6 system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
}
main { max-width: 46rem; margin: 0 auto; }
h1 { font-size: 1.9rem; line-height: 1.25; margin: 0.4em 0 0.5em; }
h2 { font-size: 1.4rem; line-height: 1.3; margin: 2.2em 0 0.6em; padding-top: 0.6em; border-top: 1px solid var(--rule); }
h3 { font-size: 1.12rem; margin: 1.6em 0 0.4em; }
p, li { overflow-wrap: anywhere; }
a { color: var(--accent); }
code { font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 0.88em;
       background: var(--code-bg); padding: 0.08em 0.3em; border-radius: 4px; }
hr { border: 0; border-top: 1px solid var(--rule); margin: 2em 0; }
blockquote { margin: 1.2em 0; padding: 0.8em 1em; border-left: 4px solid var(--rule); background: var(--box-bg); border-radius: 6px; }
blockquote.flag { border-left-color: var(--flag-edge); background: var(--flag-bg); }
blockquote.glossary { font-size: 0.95em; }
blockquote.glossary ul { padding-left: 1.1em; }
blockquote p:first-child { margin-top: 0; } blockquote p:last-child { margin-bottom: 0; }
figure { margin: 1.6em 0; }
figure img { display: block; width: 100%; height: auto; background: #ffffff; border-radius: 6px;
             padding: 6px; border: 1px solid var(--rule); }
figcaption { font-size: 0.92em; color: var(--ink2); margin-top: 0.6em; }
.table-wrap { overflow-x: auto; margin: 1.2em 0; -webkit-overflow-scrolling: touch; }
table { border-collapse: collapse; width: 100%; font-size: 0.9em; }
th, td { border-bottom: 1px solid var(--rule); padding: 0.45em 0.5em; text-align: left; vertical-align: top; }
th { background: var(--box-bg); }
ol.refs, section.refs ol { font-size: 0.88em; }
.byline { color: var(--ink2); }
"""


def inline(text):
    """Inline Markdown: code spans, bold, italic. Escapes HTML first."""
    parts = re.split(r"(`[^`]+`)", text)
    out = []
    for part in parts:
        if part.startswith("`") and part.endswith("`") and len(part) > 1:
            out.append("<code>" + html.escape(part[1:-1]) + "</code>")
            continue
        t = html.escape(part, quote=False)
        t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
        t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)
        out.append(t)
    return "".join(out)


def img_tag(alt, src):
    path = os.path.join(HERE, src)
    with open(path, "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return f'<img alt="{html.escape(alt)}" src="data:image/png;base64,{data}">'


def convert(lines):
    out, i, n = [], 0, len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()
        if not s:
            i += 1
            continue
        if s.startswith(">"):
            block = []
            while i < n and lines[i].strip().startswith(">"):
                block.append(re.sub(r"^\s*> ?", "", lines[i]))
                i += 1
            first = next((b for b in block if b.strip()), "")
            cls = "flag" if "AI-agent exercise" in first else ("glossary" if "Glossary" in first else "")
            out.append(f'<blockquote class="{cls}">' + convert(block) + "</blockquote>")
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            lvl = len(m.group(1))
            text = m.group(2)
            hid = re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")
            out.append(f'<h{lvl} id="{hid}">{inline(text)}</h{lvl}>')
            i += 1
            continue
        if re.match(r"^-{3,}$", s):
            out.append("<hr>")
            i += 1
            continue
        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)$", s)
        if m:
            alt, src = m.groups()
            i += 1
            while i < n and not lines[i].strip():
                i += 1
            cap = ""
            if i < n and lines[i].strip().startswith("**Figure"):
                cap = "<figcaption>" + inline(lines[i].strip()) + "</figcaption>"
                i += 1
            out.append("<figure>" + img_tag(alt, src) + cap + "</figure>")
            continue
        if s.startswith("|"):
            rows = []
            while i < n and lines[i].strip().startswith("|"):
                rows.append(lines[i].strip())
                i += 1
            cells = [[c.strip() for c in r.strip("|").split("|")] for r in rows]
            head, body = cells[0], [c for c in cells[1:] if not all(re.match(r"^:?-+:?$", x) for x in c)]
            t = ["<div class=\"table-wrap\"><table><thead><tr>"]
            t += [f"<th>{inline(c)}</th>" for c in head]
            t.append("</tr></thead><tbody>")
            for r in body:
                t.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>")
            t.append("</tbody></table></div>")
            out.append("".join(t))
            continue
        if re.match(r"^(\d+\.|-)\s+", s):
            ordered = bool(re.match(r"^\d+\.", s))
            items = []
            while i < n and lines[i].strip():
                ls = lines[i].strip()
                mm = re.match(r"^(\d+\.|-)\s+(.*)$", ls)
                if mm:
                    items.append(mm.group(2))
                else:
                    items[-1] += " " + ls
                i += 1
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join(f"<li>{inline(x)}</li>" for x in items) + f"</{tag}>")
            continue
        para = []
        while i < n and lines[i].strip() and not re.match(r"^(#|>|\||!\[|-{3,}$|\d+\.\s|-\s)", lines[i].strip()):
            para.append(lines[i].strip())
            i += 1
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
        else:
            i += 1
    return "\n".join(out)


def main():
    with open(SRC, encoding="utf-8") as fh:
        md = fh.read()
    lines = md.splitlines()
    title = re.sub(r"^#\s+", "", lines[0]).strip()
    body = convert(lines)
    doc = (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        f"<title>Kinase–substrate prediction briefing</title>\n"
        f"<meta name=\"description\" content=\"{html.escape(title)}\">\n"
        f"<style>{CSS}</style>\n</head>\n<body>\n<main>\n{body}\n</main>\n</body>\n</html>\n"
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"wrote {os.path.relpath(OUT, os.path.dirname(HERE))} ({len(doc) / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
