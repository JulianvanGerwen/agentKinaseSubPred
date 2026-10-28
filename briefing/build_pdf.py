"""Build briefing.pdf from briefing.md, reusing the converter in build_html.py.

Needs `xhtml2pdf` (pure Python; not part of the project environment):
    uv venv /tmp/pdfenv && uv pip install --python /tmp/pdfenv/bin/python xhtml2pdf
    /tmp/pdfenv/bin/python -I briefing/build_pdf.py
It uses a plain print stylesheet (light theme, A4), so it looks simpler than briefing.html.
"""
import html
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_html as bh  # noqa: E402
from xhtml2pdf import pisa  # noqa: E402

OUT = os.path.join(HERE, "briefing.pdf")

FONTS = os.environ.get("KSP_FONT_DIR")  # directory with DejaVuSans*.ttf (e.g. matplotlib's mpl-data/fonts/ttf)

CSS = """
@font-face { font-family: dv; src: url(@F@/DejaVuSans.ttf); }
@font-face { font-family: dv; font-weight: bold; src: url(@F@/DejaVuSans-Bold.ttf); }
@font-face { font-family: dv; font-style: italic; src: url(@F@/DejaVuSans-Oblique.ttf); }
@font-face { font-family: dv; font-weight: bold; font-style: italic; src: url(@F@/DejaVuSans-BoldOblique.ttf); }
@font-face { font-family: dvm; src: url(@F@/DejaVuSansMono.ttf); }
@page { size: A4; margin: 18mm 16mm 20mm 16mm;
        @frame footer { -pdf-frame-content: foot; bottom: 8mm; left: 16mm; right: 16mm; height: 8mm; } }
body { font-family: dv; font-size: 9pt; line-height: 1.38; color: #1a1a19; }
h1 { font-size: 19pt; line-height: 1.2; margin: 0 0 8pt 0; }
h2 { font-size: 14pt; margin: 16pt 0 4pt 0; border-top: 0.6pt solid #bbb; padding-top: 6pt; -pdf-keep-with-next: true; }
h3 { font-size: 11pt; margin: 10pt 0 3pt 0; -pdf-keep-with-next: true; }
p { margin: 0 0 6pt 0; }
li { margin: 0 0 3pt 0; }
code { font-family: dvm; font-size: 8pt; background-color: #efeeea; }
blockquote { margin: 6pt 0; padding: 5pt 8pt; background-color: #f4f3f0; border-left: 3pt solid #bbb; }
blockquote.flag { background-color: #fdf0ea; border-left: 3pt solid #eb6834; }
blockquote.glossary { font-size: 9pt; }
figure { margin: 8pt 0; }
img { width: 480pt; }
p.cap { font-size: 8.6pt; color: #52514e; margin-top: 3pt; }
table { width: 100%; font-size: 8.4pt; border-collapse: collapse; margin: 6pt 0; }
th { background-color: #f4f3f0; text-align: left; padding: 3pt; border-bottom: 0.6pt solid #999; }
td { padding: 3pt; border-bottom: 0.4pt solid #ddd; vertical-align: top; }
"""


def img_file(alt, src):
    return f'<img alt="{html.escape(alt)}" src="{os.path.join(HERE, src)}">'


def main():
    bh.img_tag = img_file
    with open(bh.SRC, encoding="utf-8") as fh:
        lines = fh.read().splitlines()
    body = bh.convert(lines)
    body = body.replace('<div class="table-wrap">', "<div>")
    body = re.sub(r"<figcaption>(.*?)</figcaption>", r'<p class="cap">\1</p>', body)
    body = body.replace("<figure>", "<div>").replace("</figure>", "</div>")
    doc = (f"<html><head><meta charset=\"utf-8\"><style>{CSS.replace("@F@", FONTS)}</style></head><body>{body}"
           '<div id="foot" style="text-align:center;font-size:8pt;color:#777">'
           "Page <pdf:pagenumber/> of <pdf:pagecount/></div></body></html>")
    with open(OUT, "wb") as fh:
        res = pisa.CreatePDF(doc, dest=fh, encoding="utf-8")
    print("errors:", res.err, "->", os.path.relpath(OUT))


if __name__ == "__main__":
    main()
