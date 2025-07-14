"""
Test calculation functions for peptide calculator.
"""

import sys
import os
import pytest
from unittest.mock import patch, MagicMock

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from src.peptide_calculator import calculate_peptide_mz, validate_openms_sequence


class TestCalculationFunctions:
    """Test cases for m/z calculation functions."""

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_basic(self, mock_charge, mock_parse, mock_pyopenms):
        mock_charge.return_value = ("PEPTIDE", 1)
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)

        result = calculate_peptide_mz("PEPTIDE", 2)

        assert "mz_ratio" in result
        assert "monoisotopic_mass" in result
        assert "molecular_formula" in result
        assert "original_sequence" in result
        assert "modified_sequence" in result
        assert "charge_state" in result
        assert "charge_source" in result
        assert "modification" in result
        assert "sequence_length" in result
        assert "aa_composition" in result
        assert "success" in result

        assert result["success"] is True
        assert result["original_sequence"] == "PEPTIDE"
        assert result["charge_state"] == 2
        assert result["charge_source"] == "From input parameter"
        assert result["sequence_length"] == 7

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_with_charge_notation(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        """Test charge override when sequence contains charge notation."""
        mock_charge.return_value = ("PEPTIDE", 3)
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 3)

        result = calculate_peptide_mz("PEPTIDE/3", 2)

        assert result["success"] is True
        assert result["charge_state"] == 3  # Should override input charge
        assert result["charge_source"] == "From sequence notation"

    @patch("src.peptide_calculator.apply_modification")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_with_dropdown_modification(
        self, mock_charge, mock_parse, mock_apply, mock_pyopenms
    ):
        """Test modification application from dropdown when ProForma parsing fails."""
        mock_charge.return_value = ("MPEPTIDE", 1)
        mock_parse.return_value = ("MPEPTIDE", "MPEPTIDE", 1)
        mock_apply.return_value = "M(Oxidation)PEPTIDE"

        # Make ProForma direct parsing fail so it goes to dropdown modification path
        mock_aa_seq_final = mock_pyopenms.AASequence.fromString.return_value
        mock_pyopenms.AASequence.fromString.side_effect = [
            Exception("Parse failed"),
            mock_aa_seq_final,
        ]

        result = calculate_peptide_mz("MPEPTIDE", 2, "Oxidation (M)")

        assert result["success"] is True
        assert result["modification"] == "Oxidation (M)"
        assert result["modified_sequence"] == "M(Oxidation)PEPTIDE"

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_with_sequence_modifications(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        """Test modification conversion from square bracket to parentheses notation."""

        def mock_from_string_side_effect(seq):
            if seq == "M[Oxidation]PEPTIDE":
                raise Exception("Parse failed")
            else:
                return mock_pyopenms.AASequence.fromString.return_value

        mock_pyopenms.AASequence.fromString.side_effect = mock_from_string_side_effect

        mock_charge.return_value = ("M[Oxidation]PEPTIDE", 1)
        mock_parse.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE", 1)

        result = calculate_peptide_mz("M[Oxidation]PEPTIDE", 2)

        assert result["success"] is True
        assert result["original_sequence"] == "MPEPTIDE"
        assert result["modified_sequence"] == "M(Oxidation)PEPTIDE"
        assert result["modification"] == "From sequence notation (converted)"

    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_proforma_sequence(self, mock_charge, mock_pyopenms):
        """Test ProForma sequence direct parsing (mass deltas)."""
        mock_charge.return_value = ("PEPTIDE[+42.0106]", 1)

        result = calculate_peptide_mz("PEPTIDE[+42.0106]", 2)

        assert result["success"] is True
        assert (
            result["modification"] == "ProForma arbitrary mass deltas (direct parsing)"
        )

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_amino_acid_composition(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        mock_charge.return_value = ("PEPTIDE", 1)
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)

        result = calculate_peptide_mz("PEPTIDE", 2)

        assert "aa_composition" in result
        aa_comp = result["aa_composition"]
        assert aa_comp["P"] == 2
        assert aa_comp["E"] == 2
        assert aa_comp["T"] == 1
        assert aa_comp["I"] == 1
        assert aa_comp["D"] == 1

    def test_calculate_peptide_mz_empty_sequence(self):
        try:
            calculate_peptide_mz("", 2)
            assert False, "Should raise ValueError for empty sequence"
        except ValueError as e:
            assert "Peptide sequence cannot be empty" in str(e)

        try:
            calculate_peptide_mz("   ", 2)
            assert False, "Should raise ValueError for whitespace-only sequence"
        except ValueError as e:
            assert "Peptide sequence cannot be empty" in str(e)

    def test_calculate_peptide_mz_invalid_charge(self):
        try:
            calculate_peptide_mz("PEPTIDE", 0)
            assert False, "Should raise ValueError for zero charge"
        except ValueError as e:
            assert "Charge state must be a positive integer" in str(e)

        try:
            calculate_peptide_mz("PEPTIDE", -1)
            assert False, "Should raise ValueError for negative charge"
        except ValueError as e:
            assert "Charge state must be a positive integer" in str(e)

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_invalid_amino_acids(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        mock_charge.return_value = ("PEPTIDEZ", 1)
        mock_parse.return_value = ("PEPTIDEZ", "PEPTIDEZ", 1)

        with pytest.raises(ValueError, match="Invalid amino acid"):
            calculate_peptide_mz("PEPTIDEZ", 2)

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_invalid_characters(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        mock_charge.return_value = ("PEPTIDE123", 1)
        mock_parse.return_value = ("PEPTIDE123", "PEPTIDE123", 1)

        with pytest.raises(ValueError, match="Invalid character"):
            calculate_peptide_mz("PEPTIDE123", 2)

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_openms_parse_failure(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        """Test failure when final OpenMS parsing fails."""
        mock_pyopenms.AASequence.fromString.side_effect = ValueError(
            "Invalid OpenMS sequence"
        )
        mock_charge.return_value = ("PEPTIDE", 1)
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)

        with pytest.raises(ValueError, match="Failed to parse modified sequence"):
            calculate_peptide_mz("PEPTIDE", 2)

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    @patch("src.peptide_calculator.parse_charge_notation")
    def test_calculate_peptide_mz_complex_sequence(
        self, mock_charge, mock_parse, mock_pyopenms
    ):
        default_return_value = mock_pyopenms.AASequence.fromString.return_value

        def from_string_side_effect(seq):
            if seq == ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]":
                raise Exception("not proforma")
            return default_return_value

        mock_pyopenms.AASequence.fromString.side_effect = from_string_side_effect
        mock_charge.return_value = (
            ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]",
            2,
        )
        mock_parse.return_value = (
            "MPEPTIDEC",
            ".(Acetyl)M(Oxidation)PEPTIDEC(Carbamidomethyl)",
            2,
        )

        result = calculate_peptide_mz(
            ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]/2", 2
        )

        assert result["success"] is True
        assert result["original_sequence"] == "MPEPTIDEC"
        assert (
            result["modified_sequence"]
            == ".(Acetyl)M(Oxidation)PEPTIDEC(Carbamidomethyl)"
        )
        assert result["charge_state"] == 2

    def test_calculate_peptide_mz_realistic_values(self, mock_pyopenms):
        """Test with realistic peptide values to ensure correct property extraction."""
        # Setup mock AASequence with specific return values
        mock_aa_seq = mock_pyopenms.AASequence.fromString.return_value
        mock_aa_seq.getMZ.return_value = 456.78
        mock_aa_seq.getMonoWeight.return_value = 911.55
        mock_aa_seq.getFormula.return_value.toString.return_value = "C40H65N11O11S"

        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_sequence_with_mods_and_charge"
            ) as mock_parse:
                mock_charge.return_value = ("MYPEPTIDE", 2)
                mock_parse.return_value = ("MYPEPTIDE", "MYPEPTIDE", 2)

                result = calculate_peptide_mz("MYPEPTIDE/2", 2)

                assert result["mz_ratio"] == 456.78
                assert result["monoisotopic_mass"] == 911.55
                assert result["molecular_formula"] == "C40H65N11O11S"

    def test_calculate_peptide_mz_x_mass_delta_fix(self, mock_pyopenms):
        """Test calculation with X[+mass] sequences - regression test for issue fix."""

        # Setup mock AASequence to fail on original sequence but succeed on corrected one
        def mock_from_string(sequence):
            if sequence == "RTAAX[+367.0537]WT":
                # This should fail, forcing the function to use our parsing logic
                raise RuntimeError(
                    "Using a mass difference to specify a modification on a residue of unknown mass is not supported"
                )
            else:
                # Return mock for the corrected sequence
                mock_aa_seq = MagicMock()
                mock_aa_seq.getMZ.return_value = 536.7144
                mock_aa_seq.getMonoWeight.return_value = 1071.4143
                mock_aa_seq.getFormula.return_value.toString.return_value = (
                    "C31H48N10O9"
                )
                return mock_aa_seq

        mock_pyopenms.AASequence.fromString.side_effect = mock_from_string

        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            mock_charge.return_value = ("RTAAX[+367.0537]WT", 1)

            result = calculate_peptide_mz("RTAAX[+367.0537]WT", 2)

            # Verify the calculation
            assert result["success"] is True
            assert result["mz_ratio"] == 536.7144
            assert result["monoisotopic_mass"] == 1071.4143
            assert result["original_sequence"] == "RTAAXWT"
            assert result["modified_sequence"] == "RTAAX[367.0537]WT"
            assert result["charge_state"] == 2
            assert result["sequence_length"] == 7

            # Verify the corrected sequence was used for PyOpenMS call
            assert mock_pyopenms.AASequence.fromString.call_count == 2


class TestValidateOpenMSSequence:
    """Test cases for OpenMS sequence validation."""

    def test_validate_openms_sequence_valid(self, mock_pyopenms):
        result = validate_openms_sequence("PEPTIDE")
        assert result is True

        result = validate_openms_sequence("M(Oxidation)PEPTIDE")
        assert result is True

    def test_validate_openms_sequence_invalid(self, mock_pyopenms):
        mock_pyopenms.AASequence.fromString.side_effect = Exception("Invalid sequence")

        result = validate_openms_sequence("INVALID")
        assert result is False

        result = validate_openms_sequence("PEPTIDE[UnknownMod]")
        assert result is False

    def test_validate_openms_sequence_empty(self, mock_pyopenms):
        """Test edge case - empty sequences might be valid in OpenMS context."""
        result = validate_openms_sequence("")
        assert result is True  # Empty might be valid for OpenMS

    def test_validate_openms_sequence_exception_handling(self, mock_pyopenms):
        mock_pyopenms.AASequence.fromString.side_effect = RuntimeError("OpenMS error")

        result = validate_openms_sequence("SEQUENCE")
        assert result is False


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
