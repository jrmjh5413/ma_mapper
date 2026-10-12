"""Figures land where they were asked to, and the table helpers filter what they say.

The save paths are the part of `plots` that changed in 0.9.0 — three functions
raised `NameError` when asked to save, the rest wrote into the installed
package directory, and `plot_heatmap` ignored a path entirely — so they are
what these tests pin down. They do not check what the figures look like.
"""

import matplotlib
import numpy as np
import pandas as pd
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ma_mapper import mapper, plots, utility  # noqa: E402


@pytest.fixture
def overlay(alignment_fasta, position_bigwig):
    """A small real overlay: the toy family against the position-valued track."""
    matrix, _, alignment = mapper.map_and_overlay(
        alignment=alignment_fasta,
        data_file=position_bigwig,
        data_format="bigwig",
        filter=False,
        return_details=True,
    )
    return matrix, alignment


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


class TestSavePaths:
    def test_heatmap_honours_a_path(self, overlay, tmp_path):
        """plot_heatmap tested `save_to_file == True`, so a path saved nothing."""
        out = tmp_path / "heatmap.png"
        plots.plot_heatmap(data=overlay[0], save_to_file=str(out))
        assert out.stat().st_size > 1000

    def test_bar_honours_a_path(self, overlay, tmp_path):
        out = tmp_path / "bar.png"
        profile = np.nan_to_num(np.nanmean(overlay[0], axis=0))
        plots.plot_bar(data=profile, save_to_file=str(out))
        assert out.stat().st_size > 1000

    def test_basecount_honours_a_path(self, overlay, tmp_path):
        """plot_basecount read an output_filepath that was never assigned."""
        out = tmp_path / "basecount.png"
        plots.plot_basecount(alignment=overlay[1], save_to_file=str(out))
        assert out.stat().st_size > 1000

    def test_colorbar_honours_a_path(self, tmp_path):
        """plot_colorbar had the same defect."""
        out = tmp_path / "colorbar.png"
        plots.plot_colorbar(cmap="Blues", vmin=0, vmax=1, save_to_file=str(out))
        assert out.stat().st_size > 500

    def test_nothing_is_written_into_the_package_directory(self, overlay, tmp_path, monkeypatch):
        """save_to_file=True used to write into the installed package directory."""
        monkeypatch.chdir(tmp_path)
        plots.plot_heatmap(data=overlay[0], save_to_file=True)
        assert (tmp_path / "overlay.png").exists()
        package_dir = __import__("ma_mapper").__path__[0]
        assert not any(name.endswith(".png") for name in __import__("os").listdir(package_dir))


class TestDrawingIntoAnAxes:
    def test_heatmap_uses_the_axes_it_is_given(self, overlay):
        fig, ax = plt.subplots()
        plots.plot_heatmap(data=overlay[0], matplot_axes=ax)
        assert len(ax.images) == 1

    def test_basecount_uses_the_axes_it_is_given(self, overlay):
        fig, ax = plt.subplots()
        plots.plot_basecount(alignment=overlay[1], matplot_axes=ax)
        assert len(ax.containers) >= 1

    def test_heatmap_rejects_nonsense_data(self):
        with pytest.raises(ValueError, match="Invalid data"):
            plots.plot_heatmap(data="not an array")


class TestRepeatmaskerPrep:
    def _table(self):
        return pd.DataFrame(
            {
                "genoName": ["chr1", "chr2", "chrUn_gl000220", "chr1", "chr19"],
                "genoStart": [10, 20, 30, 40, 50],
                "genoEnd": [100, 200, 300, 400, 500],
                "repName": ["MER11A", "L1HS", "AluY", "(CA)n", "MER11B"],
                "repClass": ["LTR", "LINE/L1", "SINE/Alu", "Simple_repeat", "LTR"],
            }
        )

    def test_drops_scaffolds_and_simple_repeats(self):
        kept = utility.repeatmasker_prep(self._table())
        assert list(kept.repName) == ["MER11A", "L1HS", "MER11B"]

    def test_another_assembly_needs_its_own_chromosome_list(self):
        """The default list is human; mouse chr-names overlap it only by accident."""
        mouse = self._table().assign(genoName=["chr1", "chr2", "GL456210.1", "chr1", "chr19"])
        default = utility.repeatmasker_prep(mouse)
        explicit = utility.repeatmasker_prep(
            mouse, main_chr=["chr1", "chr2", "chr19", "GL456210.1"]
        )
        assert len(explicit) > len(default)

    def test_an_empty_chromosome_list_keeps_everything(self):
        kept = utility.repeatmasker_prep(self._table(), main_chr=[], discard_class=[])
        assert len(kept) == 5


class TestDivergenceTablePrep:
    def test_the_filter_is_applied(self, tmp_path):
        """The filter used to be computed and then thrown away."""
        path = tmp_path / "species.tsv"
        pd.DataFrame(
            {
                "Species": ["Homo sapiens", "Mus musculus", "Absent species"],
                "Source": ["1. ok", "1. ok", "3. Zoonomia (not in alignment)"],
            }
        ).to_csv(path, sep="\t", index=False)
        kept = utility.divergence_table_prep(str(path))
        assert list(kept.Species) == ["Homo sapiens", "Mus musculus"]
        assert list(kept.meta_name) == ["Homo_sapiens", "Mus_musculus"]
        assert len(utility.divergence_table_prep(str(path), drop_unaligned=False)) == 3


class TestAgeTables:
    def test_the_age_tables_line_up(self):
        reference = utility.age_reference_table()
        lengths = {len(values) for values in reference.values()}
        assert len(lengths) == 1
        assert len(utility.age_canon()) == len(reference["age"])
        reorganised = utility.age_reference_table_reorg()
        assert len({len(values) for values in reorganised.values()}) == 1


def test_the_legacy_modules_still_import():
    """custom_cmap and gzmaf are untouched upstream code; they must at least load."""
    from ma_mapper import custom_cmap, gzmaf

    assert custom_cmap is not None and gzmaf is not None
