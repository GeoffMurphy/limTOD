"""Markdown -> PDF with typeset maths.

No LaTeX on this machine, so each $...$ / $$...$$ expression is rendered to an
SVG through matplotlib's mathtext and inlined as an <img>. Everything else goes
through python-markdown, then WeasyPrint.

Usage:
    md2pdf.py <input.md> <output.pdf> [--title "..."] [--subtitle "..."]

Needs weasyprint + markdown + matplotlib, which are NOT in the limTOD venv
(weasyprint pulls in pango/cairo bindings and is only needed for docs):

    python -m venv /tmp/docvenv && /tmp/docvenv/bin/pip install weasyprint markdown matplotlib
    /tmp/docvenv/bin/python md2pdf.py ANALYSIS_LOG.md ANALYSIS_LOG.pdf \
        --subtitle "Why each step was taken, and the maths behind it"

Gotcha worth knowing: mathtext signals a parse failure by rendering the string
LITERALLY rather than raising, so a try/except around the draw catches nothing
and the PDF silently contains raw LaTeX. That is why parse() is called
explicitly before rendering. Two things it will not accept: a bare space inside
\mathrm{}, and an expression spanning more than one line.
"""
import base64
import io
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import markdown

matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

# mathtext understands a large LaTeX subset but not these; map them over.
def _upright(m):
    """\\text{a b} / \\operatorname{a b} -> \\mathrm{a\\ b}.

    mathtext has no \\text, and a bare space inside \\mathrm{} is a parse
    error -- which it reports by silently rendering the string literally, so
    this has to be got right rather than caught.
    """
    return r"\mathrm{" + m.group(1).replace(" ", r"\ ") + "}"


FIXUPS = [
    (r"\\operatorname\{([^}]*)\}", _upright),
    (r"\\text\{([^}]*)\}", _upright),
    (r"\\!", ""),
    (r"\\,", r"\\ "),
    (r"\\;", r"\\ "),
    (r"\\quad", r"\\ \\ "),
]

_cache = {}


def render_math(expr, display, fontsize=13):
    """Render one expression to a data: URI SVG. Returns (uri, depth_ratio)."""
    key = (expr, display, fontsize)
    if key in _cache:
        return _cache[key]

    # mathtext cannot span lines, and display equations are routinely wrapped
    # in the source markdown -- collapse to one line before anything else.
    src = " ".join(expr.split())
    for pat, rep in FIXUPS:
        src = re.sub(pat, rep, src)

    # Validate first. mathtext reports a parse failure by rendering the string
    # literally rather than raising, so checking the parser explicitly is the
    # only way to notice.
    from matplotlib import mathtext
    try:
        mathtext.MathTextParser("agg").parse(f"${src}$", dpi=100, prop=None)
    except Exception as exc:
        print(f"  ! MATH PARSE FAILED, left as text: {expr[:70]}\n      {exc}")
        _cache[key] = (None, 0.0)
        return _cache[key]

    fig = plt.figure(figsize=(0.01, 0.01))
    size = fontsize * (1.25 if display else 1.0)
    t = fig.text(0, 0, f"${src}$", fontsize=size)
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())

    dpi = fig.dpi
    fig.set_size_inches(bb.width / dpi + 0.02, bb.height / dpi + 0.02)
    buf = io.StringIO()
    fig.savefig(buf, format="svg", transparent=True, bbox_inches="tight", pad_inches=0.01)
    plt.close(fig)

    uri = "data:image/svg+xml;base64," + base64.b64encode(
        buf.getvalue().encode()).decode()
    _cache[key] = (uri, bb.height / dpi)
    return _cache[key]


def substitute_math(text):
    """Replace $$..$$ then $..$ with <img> tags, skipping fenced code blocks."""
    out, in_fence = [], False
    for line in text.split("\n"):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            out.append(line)
            continue
        if in_fence:
            out.append(line)
            continue
        out.append(line)
    text = "\n".join(out)

    def disp(m):
        uri, _ = render_math(m.group(1).strip(), True)
        if uri is None:
            return m.group(0)
        return f'<div class="math-display"><img src="{uri}" alt="equation"></div>'

    def inline(m):
        uri, h = render_math(m.group(1).strip(), False)
        if uri is None:
            return m.group(0)
        return (f'<img class="math-inline" src="{uri}" alt="math" '
                f'style="height:{h*72*1.02:.1f}pt">')

    # Display first so $$ is not eaten by the single-$ pattern.
    text = re.sub(r"\$\$(.+?)\$\$", disp, text, flags=re.S)
    text = re.sub(r"(?<!\$)\$([^$\n]+?)\$(?!\$)", inline, text)
    return text


CSS = """
@page {
  size: A4; margin: 20mm 18mm 18mm 18mm;
  @bottom-center { content: counter(page) " / " counter(pages);
                   font: 8pt "DejaVu Sans"; color: #888; }
}
body { font: 10pt/1.5 "DejaVu Sans", sans-serif; color: #1a1a1a; }
h1 { font-size: 19pt; margin: 0 0 2pt; color: #10243e; line-height: 1.25; }
h2 { font-size: 13pt; margin: 20pt 0 6pt; color: #10243e;
     border-bottom: 1.5pt solid #2a78d6; padding-bottom: 3pt;
     break-after: avoid; }
h3 { font-size: 11pt; margin: 13pt 0 4pt; color: #1c4f8a; break-after: avoid; }
p, li { orphans: 2; widows: 2; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 8.6pt;
       background: #f2f4f7; padding: 0.5pt 2.5pt; border-radius: 2pt; }
pre { background: #f7f8fa; border-left: 2.5pt solid #c6d2e0; padding: 7pt 9pt;
      font-size: 8.2pt; line-height: 1.35; overflow-x: auto; break-inside: avoid; }
pre code { background: none; padding: 0; }
table { border-collapse: collapse; width: 100%; margin: 9pt 0; font-size: 8.8pt;
        break-inside: avoid; }
th { background: #eaf0f7; text-align: left; font-weight: 600; color: #10243e; }
th, td { border: 0.6pt solid #ccd6e2; padding: 3.5pt 5pt; vertical-align: top; }
tr:nth-child(even) td { background: #fafbfc; }
blockquote { margin: 8pt 0; padding: 4pt 10pt; border-left: 2.5pt solid #2a78d6;
             background: #f5f8fc; }
hr { border: none; border-top: 0.6pt solid #d8dee6; margin: 16pt 0; }
strong { color: #0d1b2e; }
.math-display { text-align: center; margin: 11pt 0; break-inside: avoid; }
.math-display img { max-width: 100%; }
.math-inline { vertical-align: -18%; }
.subtitle { color: #5a6672; font-size: 10pt; margin: 0 0 4pt; }
.built { color: #8a949e; font-size: 8pt; margin: 0 0 14pt; }
"""


def main():
    src, dst = sys.argv[1], sys.argv[2]
    title = subtitle = None
    if "--title" in sys.argv:
        title = sys.argv[sys.argv.index("--title") + 1]
    if "--subtitle" in sys.argv:
        subtitle = sys.argv[sys.argv.index("--subtitle") + 1]

    text = open(src).read()

    # Promote the leading "# " title into the styled header block.
    lines = text.split("\n")
    if lines and lines[0].startswith("# "):
        title = title or lines[0][2:].strip()
        text = "\n".join(lines[1:]).lstrip("\n")

    print("rendering maths ...")
    text = substitute_math(text)

    body = markdown.markdown(
        text, extensions=["tables", "fenced_code", "sane_lists", "attr_list"])

    from datetime import date
    header = f"<h1>{title}</h1>" if title else ""
    if subtitle:
        header += f'<p class="subtitle">{subtitle}</p>'
    header += f'<p class="built">Generated {date.today().isoformat()} from {src.split("/")[-1]}</p>'

    html = f"<html><head><meta charset='utf-8'></head><body>{header}{body}</body></html>"

    from weasyprint import HTML, CSS as WCSS
    HTML(string=html, base_url=".").write_pdf(dst, stylesheets=[WCSS(string=CSS)])
    print(f"wrote {dst}")


if __name__ == "__main__":
    main()
