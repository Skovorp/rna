from pathlib import Path

import pytest
from streamlit.testing.v1 import AppTest


APP_DIR = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("filename", "title", "sample_count"),
    [
        ("2_Ovary_paper_vs_reprocessed.py", "Ovary: paper vs reprocessed", 33),
        ("4_Atlas_paper_vs_reprocessed.py", "Tissue atlas: paper vs reprocessed", 122),
    ],
)
def test_comparison_page_renders_native_figures(
    monkeypatch, filename, title, sample_count
):
    monkeypatch.syspath_prepend(str(APP_DIR))
    # A page executed by itself has no multipage registry for st.page_link;
    # production does, so stub only that navigation element in this page test.
    monkeypatch.setattr("streamlit.page_link", lambda *args, **kwargs: None)
    app = AppTest.from_file(
        str(APP_DIR / "pages" / filename), default_timeout=45
    ).run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert [heading.value for heading in app.title] == [title]
    assert len(app.get("plotly_chart")) == 1
    assert len(app.metric) == 4
    rendered = " ".join(element.value for element in app.markdown)
    assert str(sample_count) in rendered
    assert "## TPM agreement" in rendered
    assert "## Sample PCA" not in rendered
    assert "## Sample identity" not in rendered
    assert [control.value for control in app.segmented_control] == ["TPM agreement"]
    assert any(
        button.label == "Download the standalone report"
        for button in app.download_button
    )

    for analysis, heading, table_count in (
        ("Zero transitions", "## Exact zero and non-zero transitions", 2),
        ("Sample PCA", "## Sample PCA", 0),
        ("Sample identity", "## Sample identity", 0),
    ):
        app.segmented_control[0].set_value(analysis).run()
        assert not app.exception, [exception.message for exception in app.exception]
        assert len(app.get("plotly_chart")) == 1
        assert len(app.dataframe) == table_count
        rendered = " ".join(element.value for element in app.markdown)
        assert heading in rendered
        if analysis == "Sample PCA":
            assert "500 highest-variance genes" in rendered
            assert "average of its variance within" in rendered
            assert "does not scale them to unit variance" in rendered
            assert "DESeq2-style" not in rendered
