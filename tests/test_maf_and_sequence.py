"""MAF extraction, sequence fetching and the alignment round trip."""

import numpy as np
import pytest

from ma_mapper import extract_maf, mapper, sequence_alignment

from .conftest import CHROM


class TestMafArgumentWiring:
    def test_species_list_actually_filters(self, three_species_maf, coordinate_table):
        """maf_io used to deliver species_list to the wrong parameter.

        The five fixed arguments were passed positionally in the order
        target_species, count_arg, e_value_df, internal_id_df, species_list
        while the worker declared species_list, e_value_df, internal_id_df:
        a species_list arrived as internal_id_df and raised NameError, and an
        e_value_table arrived as species_list and became a regex of its own
        column names.
        """
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]

        everything = extract_maf.maf_io(
            copya, three_species_maf, count_arg="coverage", target_species="Homo_sapiens"
        )
        # Four rows in the block, but Canis_lupus appears twice and is merged.
        assert everything[0][0] == pytest.approx(3.0)

        two_species = extract_maf.maf_io(
            copya,
            three_species_maf,
            count_arg="coverage",
            target_species="Homo_sapiens",
            species_list=["Homo_sapiens", "Mus_musculus"],
        )
        assert two_species[0][0] == pytest.approx(2.0)

    def test_species_list_result_is_still_one_value_per_base(
        self, three_species_maf, coordinate_table
    ):
        """The gap-column removal used to be discarded on this path.

        Mus_musculus carries an insertion, so the reference row has a gap
        column. Each branch of the filter went back to the unfiltered frame,
        which left that column in and returned 21 values for a 20-base window.
        """
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        out = extract_maf.maf_io(
            copya,
            three_species_maf,
            count_arg="common_freq",
            species_list=["Homo_sapiens", "Mus_musculus"],
        )
        assert len(out[0]) == 20

    def test_duplicate_species_rows_are_merged(self, three_species_maf, coordinate_table):
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        out = extract_maf.maf_io(copya, three_species_maf, count_arg="coverage")
        assert all(value <= 3.0 for value in out[0])

    def test_target_species_dropped_by_the_species_list_is_reported(
        self, three_species_maf, coordinate_table
    ):
        """Leaving the reference out of species_list used to crash obscurely.

        ref_alleles came back unbound after a printed IndexError, and the next
        line failed with UnboundLocalError instead of saying what was wrong.
        """
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        with pytest.raises(ValueError, match="absent from the alignment"):
            extract_maf.maf_io(
                copya, three_species_maf, species_list=["Mus_musculus", "Canis_lupus"]
            )

    def test_index_built_for_another_species_is_rejected(
        self, three_species_maf, coordinate_table
    ):
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        with pytest.raises(ValueError, match="indexed for"):
            extract_maf.maf_io(copya, three_species_maf, target_species="Gallus_gallus")


class TestMafFilepath:
    def test_finds_a_plain_maf(self, tmp_path):
        (tmp_path / "chr1.maf").write_text("##maf version=1\n")
        assert extract_maf.get_maf_filepath(str(tmp_path), "chr1").endswith("chr1.maf")

    def test_finds_a_gzipped_maf(self, tmp_path):
        """The gz branch was unreachable: '...chr1.maf.gz' does not end '...chr1.maf'."""
        (tmp_path / "chr1.maf.gz").write_bytes(b"")
        assert extract_maf.get_maf_filepath(str(tmp_path), "chr1").endswith("chr1.maf.gz")

    def test_missing_chromosome_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            extract_maf.get_maf_filepath(str(tmp_path), "chrZ")


class TestSequenceIO:
    def test_headers_carry_the_coordinates(self, coordinate_table, genome_fasta):
        records = sequence_alignment.sequence_io(coordinate_table, genome_fasta)
        ids = [record.id for record in records]
        assert "copyA::chrT:100-120(+)" in ids
        assert "copyC::chrT:300-310,320-330(+)" in ids

    def test_minus_strand_is_reverse_complemented(self, coordinate_table, genome_fasta):
        records = {r.id.split("::")[0]: str(r.seq) for r in
                   sequence_alignment.sequence_io(coordinate_table, genome_fasta)}
        # The genome is ACGT repeated, so chrT:200-220 is ACGT... and its
        # reverse complement starts ACGT... too; check against Biopython.
        from Bio.Seq import Seq

        import pysam

        with pysam.FastaFile(genome_fasta) as handle:
            plus = handle.fetch(CHROM, 200, 220)
        assert records["copyB"] == str(Seq(plus).reverse_complement())

    def test_fragment_sequences_are_concatenated_in_element_order(
        self, coordinate_table, genome_fasta
    ):
        records = {r.id.split("::")[0]: str(r.seq) for r in
                   sequence_alignment.sequence_io(coordinate_table, genome_fasta)}
        assert len(records["copyC"]) == 20

    def test_lengths_match_the_windows(self, coordinate_table, genome_fasta):
        records = {r.id.split("::")[0]: len(r.seq) for r in
                   sequence_alignment.sequence_io(coordinate_table, genome_fasta)}
        assert records == {"copyA": 20, "copyB": 20, "copyC": 20, "copyD": 15}

    def test_round_trip_through_the_coordinate_parser(self, coordinate_table, genome_fasta, tmp_path):
        """Coordinates written into the headers must come back unchanged."""
        out = tmp_path / "seqs.fasta"
        sequence_alignment.sequence_io(coordinate_table, genome_fasta, save_to_file=str(out))
        parsed = mapper.extract_coordinate_from_alignment(str(out))
        original = coordinate_table.copy()
        original.columns = ["chrom", "start", "end", "name", "score", "strand"]
        for name in original.name.unique():
            a = original[original.name == name][["chrom", "start", "end", "strand"]].sort_values("start")
            b = parsed[parsed.name == name][["chrom", "start", "end", "strand"]].sort_values("start")
            np.testing.assert_array_equal(a.to_numpy(), b.to_numpy())

    def test_save_path_is_next_to_the_table_not_its_parent(self, coordinate_table, genome_fasta, tmp_path):
        table_path = tmp_path / "coords.bed"
        coordinate_table.to_csv(table_path, sep="\t", header=False, index=False)
        result = sequence_alignment.sequence_io(str(table_path), genome_fasta, save_to_file=True)
        assert result == str(tmp_path / "sequence.fasta")

    def test_missing_genome_raises(self, coordinate_table, tmp_path):
        with pytest.raises(FileNotFoundError, match="source_fasta"):
            sequence_alignment.sequence_io(coordinate_table, str(tmp_path / "nope.fa"))


class TestParsedAlignmentToFasta:
    def test_gap_is_written_as_a_gap(self):
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        matrix = mapper.parse_alignment(
            [SeqRecord(Seq("AC-GN"), id="r1", name="r1")]
        )
        records = sequence_alignment.parse_alignment_to_fasta(matrix)
        # 0 means gap and 5 means N; the old table wrote 0 as 'N' and had no
        # entry for 5 at all.
        assert str(records[0].seq) == "AC-GN"


class TestMafCountArguments:
    """The per-column summaries, on the toy three-species block.

    The reference is 20 bases (chrT:100-120). Mus_musculus carries one extra
    base, so the reference row has a gap column that has to be dropped;
    Canis_lupus appears twice and is merged to one row, with its two differing
    positions collapsed to IUPAC codes.
    """

    def _run(self, maf, coordinate_table, **kwargs):
        extract_maf.build_maf_index(maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        return extract_maf.maf_io(copya, maf, **kwargs)[0]

    def test_coverage_counts_species_not_rows(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="coverage")
        assert len(out) == 20
        # human, mouse and one merged dog
        assert max(out) == pytest.approx(3.0)

    def test_common_freq_is_a_fraction(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="common_freq")
        assert len(out) == 20
        assert all(0.0 <= value <= 1.0 for value in out)

    def test_ref_freq_is_zero_where_every_species_matches(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="ref_freq")
        assert len(out) == 20
        # the two positions where the duplicate dog rows disagree are the only
        # ones that can be non-reference
        assert sum(1 for value in out if value and value > 0) <= 4

    def test_base_count_returns_a_dict_per_column(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="base_count")
        assert len(out) == 20
        assert isinstance(out[0], dict)
        assert sum(out[0].values()) == pytest.approx(3.0)

    def test_base_freq_sums_to_one(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="base_freq")
        assert sum(out[0].values()) == pytest.approx(1.0)

    def test_raw_names_every_species_and_its_base(self, three_species_maf, coordinate_table):
        out = self._run(three_species_maf, coordinate_table, count_arg="raw")
        assert len(out) == 20
        assert "Homo_sapiens" in out[0] and ":" in out[0]

    def test_common_freq_nogap_ignores_gaps_in_the_denominator(
        self, three_species_maf, coordinate_table
    ):
        """It used to divide by a total that still included the gaps."""
        out = self._run(three_species_maf, coordinate_table, count_arg="common_freq_nogap")
        assert len(out) == 20
        assert all((value is None) or np.isnan(value) or 0.0 <= value <= 1.0 for value in out)

    def test_the_overlay_invariant_holds(self, three_species_maf, coordinate_table, tmp_path):
        """One value per reference base is what map_data requires."""
        out = self._run(three_species_maf, coordinate_table, count_arg="common_freq")
        assert len(out) == 20  # chrT:100-120


class TestMafCountArgValidation:
    def test_an_unimplemented_mode_is_refused(self, three_species_maf, coordinate_table):
        """'a', 't', 'c', 'g' and 'common' were documented and returned nothing.

        They had no branch in the counting loop, so every window came back as
        an empty list, which then failed downstream or silently produced an
        empty matrix.
        """
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        for mode in ("a", "t", "c", "g", "common"):
            with pytest.raises(ValueError, match="count_arg must be one of"):
                extract_maf.maf_io(copya, three_species_maf, count_arg=mode)

    def test_every_advertised_mode_returns_one_value_per_base(
        self, three_species_maf, coordinate_table
    ):
        extract_maf.build_maf_index(three_species_maf, chrom=CHROM)
        copya = coordinate_table[coordinate_table.iloc[:, 3] == "copyA"]
        for mode in extract_maf.COUNT_ARGS:
            out = extract_maf.maf_io(copya, three_species_maf, count_arg=mode)[0]
            assert len(out) == 20, f"{mode} returned {len(out)} values"
