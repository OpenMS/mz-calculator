import sys
import os
from unittest.mock import patch, MagicMock

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from src.peptide_calculator import (
    parse_square_bracket_modifications,
    parse_charge_notation,
    parse_sequence_with_mods_and_charge,
    parse_proforma_sequence,
)


class TestParsingFunctions:
    """Test cases for sequence parsing functions."""

    def test_parse_charge_notation_slash_format(self):
        seq, charge = parse_charge_notation("PEPTIDE/2")
        assert seq == "PEPTIDE"
        assert charge == 2

        seq, charge = parse_charge_notation("PEPTIDE/3")
        assert seq == "PEPTIDE"
        assert charge == 3

        seq, charge = parse_charge_notation("SEQUENCE/5")
        assert seq == "SEQUENCE"
        assert charge == 5

    def test_parse_charge_notation_trailing_number(self):
        seq, charge = parse_charge_notation("PEPTIDE2")
        assert seq == "PEPTIDE"
        assert charge == 2

        seq, charge = parse_charge_notation("SEQUENCE3")
        assert seq == "SEQUENCE"
        assert charge == 3

        seq, charge = parse_charge_notation("PROTEIN4")
        assert seq == "PROTEIN"
        assert charge == 4

    def test_parse_charge_notation_no_charge(self):
        seq, charge = parse_charge_notation("PEPTIDE")
        assert seq == "PEPTIDE"
        assert charge == 1

        seq, charge = parse_charge_notation("SIMPLESEQUENCE")
        assert seq == "SIMPLESEQUENCE"
        assert charge == 1

    def test_parse_charge_notation_with_leading_dot(self):
        seq, charge = parse_charge_notation(".PEPTIDE/2")
        assert seq == ".PEPTIDE"
        assert charge == 2

        seq, charge = parse_charge_notation(".SEQUENCE3")
        assert seq == ".SEQUENCE"
        assert charge == 3

    def test_parse_charge_notation_invalid_charge_range(self):
        """Test invalid charge states (outside 1-20 range)."""
        seq, charge = parse_charge_notation("PEPTIDE/25")
        assert seq == "PEPTIDE/25"
        assert charge == 1

        seq, charge = parse_charge_notation("PEPTIDE/0")
        assert seq == "PEPTIDE/0"
        assert charge == 1

        seq, charge = parse_charge_notation("PEPTIDE30")
        assert seq == "PEPTIDE30"
        assert charge == 1

    def test_parse_charge_notation_edge_cases(self):
        # Empty sequence
        seq, charge = parse_charge_notation("")
        assert seq == ""
        assert charge == 1

        # Only charge notation
        seq, charge = parse_charge_notation("/2")
        assert seq == ""
        assert charge == 2

        # Whitespace handling
        seq, charge = parse_charge_notation("  PEPTIDE/2  ")
        assert seq == "PEPTIDE"
        assert charge == 2

    def test_parse_square_bracket_modifications_basic(self, mock_pyopenms):
        clean_seq, openms_seq = parse_square_bracket_modifications(
            "M[Oxidation]PEPTIDE"
        )
        assert clean_seq == "MPEPTIDE"
        assert len(openms_seq) > len(clean_seq)  # Should have modification info

    def test_parse_square_bracket_modifications_n_terminal(self, mock_pyopenms):
        clean_seq, openms_seq = parse_square_bracket_modifications("[Acetyl]PEPTIDE")
        assert clean_seq == "PEPTIDE"
        assert len(openms_seq) > len(clean_seq)

        clean_seq, openms_seq = parse_square_bracket_modifications(".[Acetyl]PEPTIDE")
        assert clean_seq == "PEPTIDE"

    def test_parse_square_bracket_modifications_c_terminal(self, mock_pyopenms):
        clean_seq, openms_seq = parse_square_bracket_modifications("PEPTIDE.[Amidated]")
        assert clean_seq == "PEPTIDE"
        assert len(openms_seq) > len(clean_seq)

    def test_parse_square_bracket_modifications_mass_delta(self, mock_pyopenms):
        """Test mass delta modification parsing (e.g., [+15.9949] for oxidation)."""
        clean_seq, openms_seq = parse_square_bracket_modifications("M[+15.9949]PEPTIDE")
        assert clean_seq == "MPEPTIDE"
        assert "[+15.9949]" in openms_seq or "15.9949" in openms_seq

        clean_seq, openms_seq = parse_square_bracket_modifications(
            "C[+57.021464]PEPTIDE"
        )
        assert clean_seq == "CPEPTIDE"

    def test_parse_square_bracket_modifications_x_mass_delta(self, mock_pyopenms):
        """Test X[+mass] parsing - special handling for unknown amino acids with mass deltas."""
        # Test the original problematic sequence
        clean_seq, openms_seq = parse_square_bracket_modifications("RTAAX[+367.0537]WT")
        assert clean_seq == "RTAAXWT"
        # The + prefix should be removed for X residues
        assert openms_seq == "RTAAX[367.0537]WT"

        # Test without + prefix (should remain unchanged)
        clean_seq, openms_seq = parse_square_bracket_modifications("RTAAX[367.0537]WT")
        assert clean_seq == "RTAAXWT"
        assert openms_seq == "RTAAX[367.0537]WT"

        # Test negative mass on X (should keep negative sign)
        clean_seq, openms_seq = parse_square_bracket_modifications("PEPX[-18.0]TIDE")
        assert clean_seq == "PEPXTIDE"
        assert openms_seq == "PEPX[-18.0]TIDE"

        # Test other mass formats on X
        clean_seq, openms_seq = parse_square_bracket_modifications("X[100]PEPTIDE")
        assert clean_seq == "XPEPTIDE"
        assert openms_seq == "X[100]PEPTIDE"

    def test_parse_square_bracket_modifications_unimod(self, mock_pyopenms):
        """Test UNIMOD notation parsing (e.g., UNIMOD:4 for carbamidomethyl)."""
        clean_seq, openms_seq = parse_square_bracket_modifications("C[UNIMOD:4]PEPTIDE")
        assert clean_seq == "CPEPTIDE"

        clean_seq, openms_seq = parse_square_bracket_modifications(
            "M[UNIMOD:35]PEPTIDE"
        )
        assert clean_seq == "MPEPTIDE"

    def test_parse_square_bracket_modifications_multiple(self, mock_pyopenms):
        clean_seq, openms_seq = parse_square_bracket_modifications(
            ".[Acetyl]M[Oxidation]PEPTIDE"
        )
        assert clean_seq == "MPEPTIDE"
        assert len(openms_seq) > len(clean_seq)

    def test_parse_square_bracket_modifications_no_modifications(self, mock_pyopenms):
        clean_seq, openms_seq = parse_square_bracket_modifications("PEPTIDE")
        assert clean_seq == "PEPTIDE"
        assert openms_seq == "PEPTIDE"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_parse_sequence_with_mods_and_charge_combined(self, mock_charge, mock_mods):
        mock_charge.return_value = ("M[Oxidation]PEPTIDE", 2)
        mock_mods.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE")

        clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge(
            "M[Oxidation]PEPTIDE/2"
        )
        assert clean_seq == "MPEPTIDE"
        assert charge == 2
        assert openms_seq == "M(Oxidation)PEPTIDE"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_parse_sequence_with_mods_and_charge_only_charge(
        self, mock_charge, mock_mods
    ):
        mock_charge.return_value = ("PEPTIDE", 3)
        mock_mods.return_value = ("PEPTIDE", "PEPTIDE")

        clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge("PEPTIDE/3")
        assert clean_seq == "PEPTIDE"
        assert charge == 3
        assert openms_seq == "PEPTIDE"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_parse_sequence_with_mods_and_charge_only_mods(
        self, mock_charge, mock_mods
    ):
        mock_charge.return_value = ("M[Oxidation]PEPTIDE", 1)
        mock_mods.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE")

        clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge(
            "M[Oxidation]PEPTIDE"
        )
        assert clean_seq == "MPEPTIDE"
        assert charge == 1
        assert openms_seq == "M(Oxidation)PEPTIDE"

    def test_parse_proforma_sequence_direct_parsing(self, mock_pyopenms):
        """Test ProForma sequence parsing with direct PyOpenMS support."""
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "PEPTIDE[+42.0106]"
        )
        assert clean_seq == "PEPTIDE"
        assert proforma_direct is True
        assert converted_seq == "PEPTIDE[+42.0106]"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_parse_proforma_sequence_fallback(self, mock_mods, mock_pyopenms):
        """Test ProForma sequence parsing with fallback to square bracket parser."""
        mock_pyopenms.AASequence.fromString.side_effect = Exception("Parse error")
        mock_mods.return_value = ("PEPTIDE", "PEPTIDE")

        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "PEPTIDE[Unknown]"
        )
        assert clean_seq == "PEPTIDE"
        assert proforma_direct is False

    def test_parse_proforma_sequence_with_leading_dot(self, mock_pyopenms):
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            ".PEPTIDE[+42.0106]"
        )
        assert clean_seq == "PEPTIDE"  # Leading dot should be removed
        assert proforma_direct is True
        assert converted_seq == "PEPTIDE[+42.0106]"


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
