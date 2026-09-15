"""Verify the downloadable sources independently reproduce displayed numbers."""

import gzip
import hashlib
from io import BytesIO
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pytest

from analysis.comparison_downloads import COMPARISONS, sample_matching


ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("study", COMPARISONS)
def test_archive_preserves_complete_sources_and_reproduces_comparison(study):
    spec = COMPARISONS[study]
    with ZipFile(ROOT / "app" / "assets" / spec.assets / f"{study}_comparison_data.zip") as archive:
        manifest = json.loads(archive.read("manifest.json"))
        for name, checksum in manifest["files_sha256"].items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == checksum
        for name, source in manifest["source_files"].items():
            original = (ROOT / source["path"]).read_bytes()
            assert hashlib.sha256(original).hexdigest() == source["sha256"]
            assert archive.read(name) == (gzip.decompress(original) if source["path"].endswith(".gz") else original)
        assert archive.read("gene_matching.tsv") == gzip.decompress(
            (ROOT / "analysis" / "results" / spec.results / "matched_genes.tsv.gz").read_bytes())

        published = pd.read_csv(archive.open("published_tpm.tsv"), sep="\t")
        reprocessed = pd.read_csv(archive.open("reprocessed_tpm.tsv"), sep="\t")
        pairs = pd.read_csv(archive.open("gene_matching.tsv"), sep="\t")
        samples = pd.read_csv(archive.open("sample_matching.tsv"), sep="\t", keep_default_na=False)
        included = samples[samples.included_in_comparison].sort_values("comparison_order")
        # Independently select source columns by their exported names, then sum
        # repeated historical rows before joining the frozen gene-pair list.
        left = published[[spec.published_id, *included.published_sample]].groupby(spec.published_id).sum()
        right = reprocessed.set_index("gene_id")[included.reprocessed_sample]
        x = np.log2(left.loc[pairs.published_id].to_numpy(float) + 1)
        y = np.log2(right.loc[pairs.reanalysis_id].to_numpy(float) + 1)
        summary = json.loads(archive.read("analysis_summary.json"))
        assert summary == json.loads((ROOT / "app" / "assets" / spec.assets / "figures.json").read_text())["summary"]
        assert np.corrcoef(x.ravel(), y.ravel())[0, 1] == pytest.approx(summary["agreement"]["pearson_log2_tpm_plus_1"], abs=1e-12)
        assert np.median(np.abs(y - x)) == pytest.approx(summary["agreement"]["median_absolute_log2_error"], abs=1e-12)
        assert np.mean(np.abs(y - x) <= 1) == pytest.approx(summary["agreement"]["fraction_abs_error_le_1"], abs=1e-12)
        assert ((x == 0) & (y > 0)).sum() == summary["zero_transitions"]["published_zero_to_reanalysis_nonzero_count"]
        assert x.shape == (summary["matched_genes"], summary["samples"])
        assert len(left) > len(pairs) and len(right) > len(pairs)
        assert not pairs.published_id.duplicated().any() and not pairs.reanalysis_id.duplicated().any()
        for name in ("used_for_published_pca", "used_for_reanalysis_pca", "used_for_joint_pca"):
            assert pairs[name].sum() == 500
        if study == "neurotranscriptome":
            assert len(published) == 16176 and len(left) == 16154
            assert set(samples.loc[~samples.included_in_comparison, "reprocessed_sample"]) == {
                "Fe_An_O_1", "Fe_Br_SF_2", "Fe_Br_SF_3"}
            assert len(reprocessed.columns) - 2 == 125
        else:
            assert samples.included_in_comparison.all() and len(samples) == 33

        # The README's copy/paste example must work using only the ZIP files.
        snippet = archive.read("README.md").decode().split("```python\n", 1)[1].split("```", 1)[0]
        def archive_open(name):
            return BytesIO(archive.read(name))
        class ArchivePandas:
            @staticmethod
            def read_csv(name, **kwargs):
                return pd.read_csv(archive_open(name), **kwargs)
        # Use the exact snippet, replacing only imports with archive-backed IO.
        snippet = snippet.replace("import pandas as pd\n", "")
        exec(snippet, {"open": archive_open, "pd": ArchivePandas, "print": lambda *args: None})


def test_sample_matching_rejects_ambiguous_reprocessed_columns():
    with pytest.raises(ValueError, match="Ambiguous sample names"):
        sample_matching("neurotranscriptome", ["Fe_An_BF_1"], ["Fe_An_BF_1", "Fe_An_BF_1.1"])
