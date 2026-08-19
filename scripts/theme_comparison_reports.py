#!/usr/bin/env python3
"""Prepare the comparison reports for embedding in the atlas.

Two transformations, both on files the generator emits standalone:

1. **Dark theme.** Rewrites the CSS block and the plotly chrome colors so an
   embedded report does not sit on a white slab inside the dark app.
   Data-carrying colors (marker/line colors of the traces) are left alone;
   only backgrounds, gridlines, axis text, and page chrome change.
2. **Slimming.** The generator inlines a full ~4.85 MB copy of plotly.js into
   each 5.2 MB report, so the actual content is under 0.4 MB. That inline copy
   is replaced with a CDN reference, cached by the browser once and reused
   across both reports and every rerun. Without this the page is unusably slow.

Idempotent: re-running on an already-prepared file is a no-op.
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPORT_DIR = ROOT / "app" / "assets" / "ovary_comparison"

# Matches app/.streamlit/config.toml.
BACKGROUND = "#0E1518"
SURFACE = "#172126"
TEXT = "#EDF5F2"
MUTED = "#9AA8A5"
BORDER = "#24323a"
GRID = "#24323a"

MARKER = "<!-- atlas-dark-theme -->"

CSS_REPLACEMENTS = [
    # Page chrome.
    ("color: #0f172a; background: #f8fafc;", f"color: {TEXT}; background: {BACKGROUND};"),
    ("background: white; border: 1px solid #e2e8f0;", f"background: {SURFACE}; border: 1px solid {BORDER};"),
    ("box-shadow: 0 1px 2px rgba(15,23,42,.04);", "box-shadow: none;"),
    # Secondary text.
    ("color: #475569;", f"color: {MUTED};"),
    ("color: #64748b;", f"color: {MUTED};"),
    # Tables.
    ("border-bottom: 1px solid #e2e8f0;", f"border-bottom: 1px solid {BORDER};"),
    ("background: #f8fafc; position: sticky;", f"background: {SURFACE}; position: sticky;"),
    ("background: #f1f5f9;", f"background: {SURFACE};"),
]

PLOTLY_REPLACEMENTS = [
    ('paper_bgcolor":"white"', f'paper_bgcolor":"{SURFACE}"'),
    ('plot_bgcolor":"white"', f'plot_bgcolor":"{SURFACE}"'),
    # Plotly's default light gridlines, in descending specificity.
    ('"gridcolor":"#C8D4E3"', f'"gridcolor":"{GRID}"'),
    ('"gridcolor":"#DFE8F3"', f'"gridcolor":"{GRID}"'),
    ('"gridcolor":"#EBF0F8"', f'"gridcolor":"{GRID}"'),
    ('"linecolor":"#EBF0F8"', f'"linecolor":"{GRID}"'),
    ('"linecolor":"#506784"', f'"linecolor":"{GRID}"'),
    ('"zerolinecolor":"#EBF0F8"', f'"zerolinecolor":"{GRID}"'),
    ('"zerolinecolor":"#C8D4E3"', f'"zerolinecolor":"{GRID}"'),
    # Plotly's default axis/title ink.
    ('"color":"#2a3f5f"', f'"color":"{TEXT}"'),
    ('"color":"#506784"', f'"color":"{MUTED}"'),
]

EXTRA_CSS = f"""
<style>{MARKER}
/* Embedded in a dark app: kill the light chrome the generator emitted. */
html, body {{ background: {BACKGROUND} !important; color: {TEXT} !important; }}
a {{ color: #53D6A5; }}
::-webkit-scrollbar {{ width: 10px; height: 10px; }}
::-webkit-scrollbar-track {{ background: {BACKGROUND}; }}
::-webkit-scrollbar-thumb {{ background: {BORDER}; border-radius: 5px; }}
</style>
"""


def externalize_plotly(html: str) -> tuple[str, float]:
    """Swap the inlined plotly.js bundle for a pinned CDN reference."""
    marker = html.find("* plotly.js")
    if marker == -1:
        return html, 0.0
    start = html.rfind("<script>", 0, marker)
    end = html.find("</script>", marker) + len("</script>")
    version_line = html[marker : html.find("\n", marker)]
    version = version_line.split("v")[-1].strip()
    saved = (end - start) / 1e6
    cdn = (
        f'<script src="https://cdn.plot.ly/plotly-{version}.min.js" '
        'charset="utf-8"></script>'
    )
    return html[:start] + cdn + html[end:], saved


def theme_report(path: Path) -> bool:
    html = path.read_text(encoding="utf-8")
    if MARKER in html:
        print(f"{path.name}: already prepared, skipping")
        return False

    before = len(html) / 1e6
    for old, new in CSS_REPLACEMENTS + PLOTLY_REPLACEMENTS:
        html = html.replace(old, new)

    html, saved = externalize_plotly(html)

    head_close = html.index("</head>")
    html = html[:head_close] + EXTRA_CSS + html[head_close:]

    path.write_text(html, encoding="utf-8")
    print(
        f"{path.name}: themed, {before:.1f} MB -> {len(html) / 1e6:.2f} MB "
        f"({saved:.1f} MB of inline plotly.js moved to CDN)"
    )
    return True


def main() -> None:
    reports = sorted(REPORT_DIR.glob("*.html"))
    if not reports:
        raise SystemExit(f"No reports found in {REPORT_DIR}")
    for report in reports:
        theme_report(report)


if __name__ == "__main__":
    main()
