"""
Test calculation functions for peptide calculator.
"""

import sys
import os
from unittest.mock import patch, MagicMock

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

# Mock pyopenms before importing to avoid dependency issues in test environment
mock_poms = MagicMock()
mock_poms.__version__ = "2.9.1"

# Mock AASequence with realistic behavior for testing
mock_aa_seq = MagicMock()
mock_aa_seq.getMZ.return_value = 500.25
mock_aa_seq.getMonoWeight.return_value = 999.49
mock_aa_seq.getFormula.return_value.toString.return_value = "C43H66N12O12"
mock_poms.AASequence.fromString.return_value = mock_aa_seq

sys.modules["pyopenms"] = mock_poms

from src.peptide_calculator import calculate_peptide_mz, validate_openms_sequence


class TestCalculationFunctions:
    """Test cases for m/z calculation functions."""

    def test_calculate_peptide_mz_basic(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
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

    def test_calculate_peptide_mz_with_charge_notation(self):
        """Test charge override when sequence contains charge notation."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
                    mock_charge.return_value = ("PEPTIDE", 3)
                    mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 3)

                    result = calculate_peptide_mz("PEPTIDE/3", 2)

                    assert result["success"] is True
                    assert result["charge_state"] == 3  # Should override input charge
                    assert result["charge_source"] == "From sequence notation"

    def test_calculate_peptide_mz_with_dropdown_modification(self):
        """Test modification application from dropdown when ProForma parsing fails."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            # Make ProForma direct parsing fail so it goes to dropdown modification path
            mock_poms_failing = MagicMock()
            mock_poms_failing.AASequence.fromString.side_effect = Exception(
                "Parse failed"
            )
            mock_get_poms.return_value = mock_poms_failing

            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
                    with patch(
                        "src.peptide_calculator.apply_modification"
                    ) as mock_apply:
                        mock_charge.return_value = ("MPEPTIDE", 1)
                        mock_parse.return_value = ("MPEPTIDE", "MPEPTIDE", 1)
                        mock_apply.return_value = "M(Oxidation)PEPTIDE"

                        # Now set up successful AASequence for the final parsing
                        mock_aa_seq_final = MagicMock()
                        mock_aa_seq_final.getMZ.return_value = 500.25
                        mock_aa_seq_final.getMonoWeight.return_value = 999.49
                        mock_aa_seq_final.getFormula.return_value.toString.return_value = (
                            "C43H66N12O12"
                        )
                        mock_poms_failing.AASequence.fromString.side_effect = [
                            Exception("Parse failed"),
                            mock_aa_seq_final,
                        ]

                        result = calculate_peptide_mz("MPEPTIDE", 2, "Oxidation (M)")

                        assert result["success"] is True
                        assert result["modification"] == "Oxidation (M)"
                        assert result["modified_sequence"] == "M(Oxidation)PEPTIDE"

    def test_calculate_peptide_mz_with_sequence_modifications(self):
        """Test modification conversion from square bracket to parentheses notation."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            # Create a local mock to control AASequence behavior for this test
            local_mock_poms = MagicMock()
            local_mock_poms.__version__ = "2.9.1"

            # Make the first call to fromString fail (for ProForma parsing)
            # but the second call succeed (for the modified sequence parsing)
            def mock_from_string_side_effect(seq):
                if seq == "M[Oxidation]PEPTIDE":
                    raise Exception("Parse failed")
                else:
                    # Return the successful mock for the modified sequence
                    return mock_aa_seq

            local_mock_poms.AASequence.fromString.side_effect = (
                mock_from_string_side_effect
            )
            mock_get_poms.return_value = local_mock_poms

            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
                    mock_charge.return_value = ("M[Oxidation]PEPTIDE", 1)
                    mock_parse.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE", 1)

                    result = calculate_peptide_mz("M[Oxidation]PEPTIDE", 2)

                    assert result["success"] is True
                    assert result["original_sequence"] == "MPEPTIDE"
                    assert result["modified_sequence"] == "M(Oxidation)PEPTIDE"
                    assert (
                        result["modification"] == "From sequence notation (converted)"
                    )

    def test_calculate_peptide_mz_proforma_sequence(self):
        """Test ProForma sequence direct parsing (mass deltas)."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            # Mock successful ProForma parsing
            mock_aa_seq_proforma = MagicMock()
            mock_aa_seq_proforma.getMZ.return_value = 500.25
            mock_aa_seq_proforma.getMonoWeight.return_value = 999.49
            mock_aa_seq_proforma.getFormula.return_value.toString.return_value = (
                "C43H66N12O12"
            )
            mock_get_poms.return_value.AASequence.fromString.return_value = (
                mock_aa_seq_proforma
            )

            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                mock_charge.return_value = ("PEPTIDE[+42.0106]", 1)

                result = calculate_peptide_mz("PEPTIDE[+42.0106]", 2)

                assert result["success"] is True
                assert (
                    result["modification"]
                    == "ProForma arbitrary mass deltas (direct parsing)"
                )

    def test_calculate_peptide_mz_amino_acid_composition(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
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

    def test_calculate_peptide_mz_invalid_amino_acids(self):
        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_sequence_with_mods_and_charge"
            ) as mock_parse:
                mock_charge.return_value = ("PEPTIDEZ", 1)
                mock_parse.return_value = ("PEPTIDEZ", "PEPTIDEZ", 1)

                try:
                    calculate_peptide_mz("PEPTIDEZ", 2)
                    assert False, "Should raise ValueError for invalid amino acids"
                except ValueError as e:
                    assert "Invalid amino acid" in str(e)
                    assert "Z" in str(e)

    def test_calculate_peptide_mz_invalid_characters(self):
        with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
            with patch(
                "src.peptide_calculator.parse_sequence_with_mods_and_charge"
            ) as mock_parse:
                mock_charge.return_value = ("PEPTIDE123", 1)
                mock_parse.return_value = ("PEPTIDE123", "PEPTIDE123", 1)

                try:
                    calculate_peptide_mz("PEPTIDE123", 2)
                    assert False, "Should raise ValueError for invalid characters"
                except ValueError as e:
                    assert "Invalid character" in str(e)

    def test_calculate_peptide_mz_openms_parse_failure(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
                    mock_charge.return_value = ("PEPTIDE", 1)
                    mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)
                    mock_get_poms.return_value.AASequence.fromString.side_effect = (
                        Exception("OpenMS parse error")
                    )

                    try:
                        calculate_peptide_mz("PEPTIDE", 2)
                        assert (
                            False
                        ), "Should raise ValueError when OpenMS parsing fails"
                    except ValueError as e:
                        assert "Failed to parse modified sequence" in str(e)

    def test_calculate_peptide_mz_complex_sequence(self):
        """Test sequence with N-terminal, modification, and C-terminal modifications."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value = mock_poms
            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
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
                        ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]/2", 1
                    )

                    assert result["success"] is True
                    assert result["charge_state"] == 2
                    assert result["charge_source"] == "From sequence notation"
                    assert len(result["modified_sequence"]) > len(
                        result["original_sequence"]
                    )

    def test_calculate_peptide_mz_realistic_values(self):
        """Test with realistic mass spec values for verification."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            # Use more realistic mock values
            mock_aa_seq_realistic = MagicMock()
            mock_aa_seq_realistic.getMZ.return_value = (
                444.2345  # Realistic m/z for PEPTIDE +2
            )
            mock_aa_seq_realistic.getMonoWeight.return_value = (
                886.4576  # Realistic mass
            )
            mock_aa_seq_realistic.getFormula.return_value.toString.return_value = (
                "C43H66N12O12"
            )
            mock_get_poms.return_value.AASequence.fromString.return_value = (
                mock_aa_seq_realistic
            )

            with patch("src.peptide_calculator.parse_charge_notation") as mock_charge:
                with patch(
                    "src.peptide_calculator.parse_sequence_with_mods_and_charge"
                ) as mock_parse:
                    mock_charge.return_value = ("PEPTIDE", 1)
                    mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)

                    result = calculate_peptide_mz("PEPTIDE", 2)

                    assert result["mz_ratio"] == 444.2345
                    assert result["monoisotopic_mass"] == 886.4576
                    assert result["molecular_formula"] == "C43H66N12O12"


class TestValidateOpenMSSequence:
    """Test cases for OpenMS sequence validation."""

    def test_validate_openms_sequence_valid(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value.AASequence.fromString.return_value = MagicMock()

            result = validate_openms_sequence("PEPTIDE")
            assert result is True

            result = validate_openms_sequence("M(Oxidation)PEPTIDE")
            assert result is True

    def test_validate_openms_sequence_invalid(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value.AASequence.fromString.side_effect = Exception(
                "Invalid sequence"
            )

            result = validate_openms_sequence("INVALID")
            assert result is False

            result = validate_openms_sequence("PEPTIDE[UnknownMod]")
            assert result is False

    def test_validate_openms_sequence_empty(self):
        """Test edge case - empty sequences might be valid in OpenMS context."""
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value.AASequence.fromString.return_value = MagicMock()

            result = validate_openms_sequence("")
            assert result is True  # Empty might be valid for OpenMS

    def test_validate_openms_sequence_exception_handling(self):
        with patch("src.peptide_calculator._get_pyopenms") as mock_get_poms:
            mock_get_poms.return_value.AASequence.fromString.side_effect = RuntimeError(
                "OpenMS error"
            )

            result = validate_openms_sequence("SEQUENCE")
            assert result is False


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
