"""Render D2_report_2026EEY7565.md to PDF: Markdown -> HTML (MathJax) -> headless Edge.

Run from this folder (Git Bash):
    PYTHONUTF8=1 uv run --no-project --with markdown-it-py --with mdit-py-plugins python build_pdf.py
"""
import html
import pathlib
import re
import subprocess

from markdown_it import MarkdownIt
from mdit_py_plugins.dollarmath import dollarmath_plugin

HERE = pathlib.Path(__file__).parent
NAME = "D2_report_2026EEY7565"
EDGE = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"

CSS = """
@page { size: A4; margin: 20mm 19mm 22mm 19mm;
        @bottom-center { content: counter(page); font: 9pt Georgia, serif; color: #666; } }
body { font: 10.5pt/1.5 Georgia, 'Times New Roman', serif; color: #1a1a1a; max-width: 172mm; margin: 0 auto; }
h1 { font-size: 20pt; color: #1F5FA8; margin: 0 0 4pt; line-height: 1.2; }
h2 { font-size: 14pt; color: #1F5FA8; margin: 20pt 0 6pt; border-bottom: 0.6pt solid #C9CED6; padding-bottom: 2pt; }
h3 { font-size: 11.5pt; margin: 14pt 0 4pt; }
h2, h3 { break-after: avoid; }
p { margin: 0 0 7pt; text-align: justify; hyphens: auto; }
img { display: block; max-width: 100%; margin: 8pt auto 3pt; }
p:has(> img) { break-inside: avoid; break-after: avoid; text-align: center; }
p:has(> img) + p { font-size: 9pt; line-height: 1.35; color: #333; margin: 0 6mm 12pt; break-inside: avoid; }
table { border-collapse: collapse; margin: 6pt auto 10pt; font-size: 9.3pt; break-inside: avoid; }
th { border-top: 1pt solid #222; border-bottom: 0.6pt solid #222; padding: 3pt 5pt; text-align: left; }
td { padding: 2.2pt 5pt; }
tr:last-child td { border-bottom: 1pt solid #222; }
blockquote { background: #F4F6F9; border: 0.5pt solid #C9CED6; border-left: 3pt solid #1F5FA8;
             margin: 8pt 0 10pt; padding: 6pt 10pt 1pt; break-inside: avoid; }
code { font: 9pt Consolas, monospace; background: #F4F6F9; padding: 0 2pt; }
pre { background: #F4F6F9; border: 0.5pt solid #C9CED6; padding: 6pt 8pt; break-inside: avoid; }
pre code { background: none; padding: 0; }
.appendix table { font-size: 7.2pt; line-height: 1.25; } .appendix td, .appendix th { padding: 0.3pt 5pt; text-align: right; }
.appendix p { break-after: avoid; break-inside: avoid; margin-top: 8pt; }
mjx-container[display="true"] { margin: 6pt 0 !important; }
"""


def render_math(content, opts):
    left, right = (r"\[", r"\]") if opts["display_mode"] else (r"\(", r"\)")
    return left + html.escape(content) + right


text = (HERE / f"{NAME}.md").read_text(encoding="utf-8")
text = re.sub(r"^> \[!\w+\] *(.*)$", r"> **\1**", text, flags=re.M)      # Obsidian callout title -> bold line
md = MarkdownIt("commonmark").enable("table").use(dollarmath_plugin, renderer=render_math)
body = md.render(text)
body = body.replace('<h2>Appendix A', '<div class="appendix"><h2>Appendix A', 1)
body = body.replace('<h2>Appendix B', '</div><h2>Appendix B', 1)

page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>{NAME}</title>
<style>{CSS}</style>
<script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
</head><body>{body}</body></html>"""
html_path = HERE / f"{NAME}.html"
html_path.write_text(page, encoding="utf-8")

subprocess.run([EDGE, "--headless=new", "--disable-gpu", "--no-pdf-header-footer",
                "--virtual-time-budget=30000", f"--print-to-pdf={HERE / (NAME + '.pdf')}",
                html_path.as_uri()], check=True)
html_path.unlink()                                   # the HTML is only an intermediate
print("wrote", HERE / f"{NAME}.pdf")
