from pathlib import Path

from streamlit.testing.v1 import AppTest


CHEATSHEET = Path(__file__).resolve().parents[1] / "pages" / "1_Mosquito_cheatsheet.py"
ANATOMY_IMAGE = Path(__file__).resolve().parents[1] / "assets" / "mosquito_anatomy.png"


def test_cheatsheet_renders_dataset_terms(monkeypatch):
    assert ANATOMY_IMAGE.exists()
    assert ANATOMY_IMAGE.stat().st_size > 100_000
    monkeypatch.syspath_prepend(str(CHEATSHEET.parents[1]))
    app = AppTest.from_file(str(CHEATSHEET), default_timeout=30).run()

    assert not app.exception, [exception.message for exception in app.exception]
    assert [title.value for title in app.title] == ["Mosquito anatomy & stages"]

    rendered = " ".join(
        str(element.value)
        for element_type in ("markdown", "caption", "info", "code")
        for element in getattr(app, element_type)
    ).casefold()
    for term in (
        "antenna",
        "compound eyes",
        "palps",
        "proboscis",
        "rostrum",
        "wings",
        "halteres",
        "forelegs",
        "femur",
        "tibia",
        "tarsus",
        "abdominal tip",
        "ovaries",
        "egg",
        "larva",
        "pupa",
        "adult",
        "blood-fed",
        "gravid",
        "eggs retained",
        "post-blood-meal",
    ):
        assert term in rendered
