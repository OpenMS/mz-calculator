"""
Test analysis and utility functions for peptide calculator.
Covers functions not yet tested in other test files.
"""

import sys
import os
from unittest.mock import patch, MagicMock

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from src.peptide_calculator import (
    analyze_peptide_sequence,
    parse_proforma_sequence,
    get_cached_modifications,
    get_cached_examples,
    SequenceAnalysis,
    ERROR_MESSAGES,
)


class TestAnalyzePeptideSequence:
    """Test cases for the unified sequence analysis function."""

    def test_analyze_empty_sequence(self):
        analysis = analyze_peptide_sequence("")

        assert analysis.is_valid is False
        assert analysis.error_message == ERROR_MESSAGES["empty_sequence"]
        assert analysis.modification == "None"
        assert analysis.charge == 2
        assert analysis.clean_sequence == ""

    def test_analyze_whitespace_sequence(self):
        analysis = analyze_peptide_sequence("   ")

        assert analysis.is_valid is False
        assert analysis.error_message == ERROR_MESSAGES["empty_sequence"]

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_basic_valid_sequence(self, mock_parse, mock_detect, mock_validate):
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 1)
        mock_detect.return_value = "None"
        mock_validate.return_value = (True, "PEPTIDE")

        analysis = analyze_peptide_sequence("PEPTIDE")

        assert analysis.is_valid is True
        assert analysis.error_message is None
        assert analysis.modification == "None"
        assert analysis.modification_detected is False
        assert analysis.charge == 2  # Default charge
        assert analysis.charge_detected is False
        assert analysis.clean_sequence == "PEPTIDE"

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_sequence_with_charge_detection(
        self, mock_parse, mock_detect, mock_validate
    ):
        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 3)
        mock_detect.return_value = "None"
        mock_validate.return_value = (True, "PEPTIDE")

        analysis = analyze_peptide_sequence("PEPTIDE/3")

        assert analysis.is_valid is True
        assert analysis.charge == 3
        assert analysis.charge_detected is True

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_sequence_with_modification_detection(
        self, mock_parse, mock_detect, mock_validate
    ):
        mock_parse.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE", 1)
        mock_detect.return_value = "Oxidation (M)"
        mock_validate.return_value = (True, "MPEPTIDE")

        analysis = analyze_peptide_sequence("M[Oxidation]PEPTIDE")

        assert analysis.is_valid is True
        assert analysis.modification == "Oxidation (M)"
        assert analysis.modification_detected is True

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_sequence_with_both_charge_and_modification(
        self, mock_parse, mock_detect, mock_validate
    ):
        mock_parse.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE", 2)
        mock_detect.return_value = "Oxidation (M)"
        mock_validate.return_value = (True, "MPEPTIDE")

        analysis = analyze_peptide_sequence("M[Oxidation]PEPTIDE/2")

        assert analysis.is_valid is True
        assert analysis.modification == "Oxidation (M)"
        assert analysis.modification_detected is True
        assert analysis.charge == 2
        assert analysis.charge_detected is True

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_invalid_sequence(self, mock_parse, mock_detect, mock_validate):
        mock_parse.return_value = ("PEPTIDEZ", "PEPTIDEZ", 1)
        mock_detect.return_value = "None"
        mock_validate.return_value = (False, "PEPTIDEZ")

        analysis = analyze_peptide_sequence("PEPTIDEZ")

        assert analysis.is_valid is False
        assert analysis.error_message == ERROR_MESSAGES["invalid_amino_acid"]

    @patch("src.peptide_calculator.validate_peptide_sequence")
    @patch("src.peptide_calculator.detect_modification_from_sequence")
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_empty_clean_sequence(self, mock_parse, mock_detect, mock_validate):
        """Test analysis when clean sequence is empty after parsing modifications."""
        mock_parse.return_value = ("", "", 1)
        mock_detect.return_value = "None"
        mock_validate.return_value = (True, "")

        analysis = analyze_peptide_sequence("[Acetyl]")

        assert analysis.is_valid is False
        assert analysis.error_message == ERROR_MESSAGES["invalid_sequence_length"]

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_analyze_sequence_exception_handling(self, mock_parse):
        mock_parse.side_effect = Exception("Parse error")

        analysis = analyze_peptide_sequence("INVALID")

        assert analysis.is_valid is False
        assert analysis.error_message is not None
        assert "unexpected error" in analysis.error_message.lower()

    def test_analyze_sequence_defaults(self):
        analysis = SequenceAnalysis()

        assert analysis.modification == "None"
        assert analysis.modification_detected is False
        assert analysis.charge == 2
        assert analysis.charge_detected is False
        assert analysis.is_valid is True
        assert analysis.clean_sequence == ""
        assert analysis.error_message is None


class TestParseProFormaSequence:
    """Test cases for ProForma sequence parsing."""

    def test_parse_proforma_direct_success(self, mock_pyopenms):
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "PEPTIDE[+42.0106]"
        )

        assert clean_seq == "PEPTIDE"
        assert converted_seq == "PEPTIDE[+42.0106]"
        assert proforma_direct is True

    def test_parse_proforma_with_leading_dot(self, mock_pyopenms):
        """Test ProForma parsing with leading dot (common in mass spec notation)."""
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            ".PEPTIDE[+42.0106]"
        )

        assert clean_seq == "PEPTIDE"
        assert converted_seq == "PEPTIDE[+42.0106]"
        assert proforma_direct is True

    def test_parse_proforma_direct_failure_fallback(self, mock_pyopenms):
        """Test fallback to traditional conversion when direct ProForma parsing fails."""
        # Make direct parsing fail
        mock_pyopenms.AASequence.fromString.side_effect = Exception("Parse error")

        with patch(
            "src.peptide_calculator.parse_square_bracket_modifications"
        ) as mock_parse:
            mock_parse.return_value = ("PEPTIDE", "PEPTIDE(Oxidation)")

            clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
                "PEPTIDE[Oxidation]"
            )

            assert clean_seq == "PEPTIDE"
            assert converted_seq == "PEPTIDE(Oxidation)"
            assert proforma_direct is False

    def test_parse_proforma_complex_mass_deltas(self, mock_pyopenms):
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "EM[+15.9949]EVEES[-79.9663]PEK"
        )

        assert clean_seq == "EMEVEESPEK"
        assert converted_seq == "EM[+15.9949]EVEES[-79.9663]PEK"
        assert proforma_direct is True

    def test_parse_proforma_scientific_notation(self, mock_pyopenms):
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "PEPTIDE[+1.5e2]"
        )

        assert clean_seq == "PEPTIDE"
        assert converted_seq == "PEPTIDE[+1.5e2]"
        assert proforma_direct is True

    def test_parse_proforma_multiple_modifications(self, mock_pyopenms):
        clean_seq, converted_seq, proforma_direct = parse_proforma_sequence(
            "K[+28.0313]PEPTIDER[-10.0086]"
        )

        assert clean_seq == "KPEPTIDER"
        assert converted_seq == "K[+28.0313]PEPTIDER[-10.0086]"
        assert proforma_direct is True


class TestCachingFunctions:
    """Test cases for caching utility functions."""

    def test_get_cached_modifications(self):
        with patch(
            "src.peptide_calculator.get_supported_modifications"
        ) as mock_get_mods:
            expected_mods = ["None", "Oxidation (M)", "Carbamidomethyl (C)"]
            mock_get_mods.return_value = expected_mods

            result = get_cached_modifications()

            assert result == expected_mods
            mock_get_mods.assert_called_once()

    def test_get_cached_examples(self):
        with patch(
            "src.peptide_calculator.get_square_bracket_examples"
        ) as mock_get_examples:
            expected_examples = {
                "M[Oxidation]PEPTIDE": "Methionine oxidation",
                "PEPTIDE/2": "Charge state 2",
            }
            mock_get_examples.return_value = expected_examples

            result = get_cached_examples()

            assert result == expected_examples
            mock_get_examples.assert_called_once()

    def test_cached_functions_return_correct_types(self):
        modifications = get_cached_modifications()
        assert isinstance(modifications, list)

        examples = get_cached_examples()
        assert isinstance(examples, dict)

    def test_cached_modifications_content(self):
        """Test content of cached modifications - verifies expected modifications are present."""
        modifications = get_cached_modifications()

        # Should contain basic modifications
        assert "None" in modifications
        assert "Oxidation (M)" in modifications
        assert "Carbamidomethyl (C)" in modifications

        assert len(modifications) > 0
        assert all(isinstance(mod, str) for mod in modifications)

    def test_cached_examples_content(self):
        """Test content of cached examples - ensures coverage of different example types."""
        examples = get_cached_examples()

        # Should contain various example types
        has_oxidation = any("Oxidation" in key for key in examples.keys())
        has_charge = any("/" in key for key in examples.keys())
        has_mass_delta = any("[+" in key for key in examples.keys())

        assert has_oxidation
        assert has_charge or has_mass_delta  # At least one of these should be present

        assert all(isinstance(desc, str) for desc in examples.values())
        assert all(isinstance(seq, str) for seq in examples.keys())
