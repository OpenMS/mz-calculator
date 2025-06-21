import sys
import os
from unittest.mock import patch, MagicMock

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

# Mock pyopenms before importing - required since pyopenms may not be available in test environment
mock_poms = MagicMock()
mock_poms.__version__ = "2.9.1"
mock_mod_db = MagicMock()
mock_mod = MagicMock()
mock_mod.getId.return_value = "Oxidation"
mock_mod_db.getModification.return_value = mock_mod
mock_poms.ModificationsDB.return_value = mock_mod_db
sys.modules["pyopenms"] = mock_poms

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

    def test_parse_square_bracket_modifications_basic(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "M[Oxidation]PEPTIDE"
                )
                assert clean_seq == "MPEPTIDE"
                assert len(openms_seq) > len(clean_seq)  # Should have modification info

    def test_parse_square_bracket_modifications_n_terminal(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "[Acetyl]PEPTIDE"
                )
                assert clean_seq == "PEPTIDE"
                assert len(openms_seq) > len(clean_seq)

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    ".[Acetyl]PEPTIDE"
                )
                assert clean_seq == "PEPTIDE"

    def test_parse_square_bracket_modifications_c_terminal(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "PEPTIDE.[Amidated]"
                )
                assert clean_seq == "PEPTIDE"
                assert len(openms_seq) > len(clean_seq)

    def test_parse_square_bracket_modifications_mass_delta(self):
        """Test mass delta modification parsing (e.g., [+15.9949] for oxidation)."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "M[+15.9949]PEPTIDE"
                )
                assert clean_seq == "MPEPTIDE"
                assert "[+15.9949]" in openms_seq or "15.9949" in openms_seq

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "C[+57.021464]PEPTIDE"
                )
                assert clean_seq == "CPEPTIDE"

    def test_parse_square_bracket_modifications_unimod(self):
        """Test UNIMOD notation parsing (e.g., UNIMOD:4 for carbamidomethyl)."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "C[UNIMOD:4]PEPTIDE"
                )
                assert clean_seq == "CPEPTIDE"

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    "M[UNIMOD:35]PEPTIDE"
                )
                assert clean_seq == "MPEPTIDE"

    def test_parse_square_bracket_modifications_multiple(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications(
                    ".[Acetyl]M[Oxidation]PEPTIDE"
                )
                assert clean_seq == "MPEPTIDE"
                assert len(openms_seq) > len(clean_seq)

    def test_parse_square_bracket_modifications_no_modifications(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator._get_pyopenms_mod_db") as mock_get_db:
                mock_get_db.return_value = mock_mod_db

                clean_seq, openms_seq = parse_square_bracket_modifications("PEPTIDE")
                assert clean_seq == "PEPTIDE"
                assert openms_seq == "PEPTIDE"

    def test_parse_sequence_with_mods_and_charge_combined(self):
        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_square_bracket_modifications"
            ) as mock_mods:
                mock_charge.return_value = ("M[Oxidation]PEPTIDE", 2)
                mock_mods.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE")

                clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge(
                    "M[Oxidation]PEPTIDE/2"
                )
                assert clean_seq == "MPEPTIDE"
                assert charge == 2
                assert openms_seq == "M(Oxidation)PEPTIDE"

    def test_parse_sequence_with_mods_and_charge_only_charge(self):
        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_square_bracket_modifications"
            ) as mock_mods:
                mock_charge.return_value = ("PEPTIDE", 3)
                mock_mods.return_value = ("PEPTIDE", "PEPTIDE")

                clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge(
                    "PEPTIDE/3"
                )
                assert clean_seq == "PEPTIDE"
                assert charge == 3
                assert openms_seq == "PEPTIDE"

    def test_parse_sequence_with_mods_and_charge_only_mods(self):
        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_square_bracket_modifications"
            ) as mock_mods:
                mock_charge.return_value = ("M[Oxidation]PEPTIDE", 1)
                mock_mods.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE")

                clean_seq, openms_seq, charge = parse_sequence_with_mods_and_charge(
                    "M[Oxidation]PEPTIDE"
                )
                assert clean_seq == "MPEPTIDE"
                assert charge == 1
                assert openms_seq == "M(Oxidation)PEPTIDE"

    def test_parse_proforma_sequence_direct_parsing(self):
        """Test ProForma sequence parsing with direct PyOpenMS support."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_aa_seq = MagicMock()
            mock_get_poms.return_value.AASequence.fromString.return_value = mock_aa_seq

            clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
                "PEPTIDE[+42.0106]"
            )
            assert clean_seq == "PEPTIDE"
            assert proforma_direct is True
            assert converted_seq == "PEPTIDE[+42.0106]"

    def test_parse_proforma_sequence_fallback(self):
        """Test ProForma sequence parsing with fallback to square bracket parser."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value.AASequence.fromString.side_effect = Exception(
                "Parse error"
            )
            with patch(
                "src.peptide_calculator.parse_square_bracket_modifications"
            ) as mock_parse:
                mock_parse.return_value = ("PEPTIDE", "PEPTIDE")

                clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
                    "PEPTIDE[Unknown]"
                )
                assert clean_seq == "PEPTIDE"
                assert proforma_direct is False

    def test_parse_proforma_sequence_with_leading_dot(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_aa_seq = MagicMock()
            mock_get_poms.return_value.AASequence.fromString.return_value = mock_aa_seq

            clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
                ".PEPTIDE[+42.0106]"
            )
            assert clean_seq == "PEPTIDE"  # Leading dot should be removed
            assert proforma_direct is True


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
