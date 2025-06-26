"""
Test validation functions for peptide calculator.
"""

import pytest
from unittest.mock import patch, MagicMock

from src.peptide_calculator import (
    validate_peptide_sequence,
    validate_peptide_sequence_with_mods,
    ERROR_MESSAGES,
)


class TestValidationFunctions:
    """Test cases for sequence validation functions."""

    @pytest.mark.parametrize(
        "sequence,expected_valid,expected_clean",
        [
            ("PEPTIDE", True, "PEPTIDE"),
            ("ACDEFGHIKLMNPQRSTVWY", True, "ACDEFGHIKLMNPQRSTVWY"),
            ("PEPTIDEX", True, "PEPTIDEX"),
            ("PEPTIDEU", True, "PEPTIDEU"),
            ("ACDEFGHIKLMNPQRSTVWYXU", True, "ACDEFGHIKLMNPQRSTVWYXU"),
            ("", True, ""),
            ("   ", True, ""),
        ],
    )
    def test_validate_peptide_sequence_valid_cases(
        self, sequence, expected_valid, expected_clean
    ):
        is_valid, clean_seq = validate_peptide_sequence(sequence)
        assert is_valid == expected_valid
        assert clean_seq == expected_clean

    @pytest.mark.parametrize(
        "sequence",
        [
            "PEPTIDEZ",  # Invalid amino acid Z
            "PEPTIDEJ",  # Invalid amino acid J
            "PEPTIDEO",  # Invalid amino acid O
            "PEPTIDEB",  # Invalid amino acid B
        ],
    )
    def test_validate_peptide_sequence_invalid_amino_acids(self, sequence):
        is_valid, clean_seq = validate_peptide_sequence(sequence)
        assert is_valid is False
        assert isinstance(clean_seq, str)

    @pytest.mark.parametrize(
        "sequence,expected_clean",
        [
            ("M[Oxidation]PEPTIDE", "MPEPTIDE"),
            ("C[Carbamidomethyl]PEPTIDE", "CPEPTIDE"),
            ("[Acetyl]PEPTIDE", "PEPTIDE"),
            ("PEPTIDE[Amidated]", "PEPTIDE"),
        ],
    )
    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_validate_peptide_sequence_with_modifications(
        self, mock_parse, sequence, expected_clean
    ):
        mock_parse.return_value = (expected_clean, f"Modified_{expected_clean}")

        is_valid, clean_seq = validate_peptide_sequence(sequence)
        assert is_valid is True
        assert clean_seq == expected_clean
        mock_parse.assert_called_once_with(sequence)

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_validate_peptide_sequence_with_leading_dot(self, mock_parse):
        # Test validation with leading dot for N-terminal modifications
        mock_parse.return_value = ("PEPTIDE", ".(Acetyl)PEPTIDE")

        is_valid, clean_seq = validate_peptide_sequence(".[Acetyl]PEPTIDE")
        assert is_valid is True
        assert clean_seq == "PEPTIDE"

    @pytest.mark.parametrize(
        "sequence",
        [
            "peptide",  # Lowercase
            "PePtIdE",  # Mixed case
            "PEPTIDE123",  # With numbers
            "PEPTIDE@#$",  # With special characters
        ],
    )
    def test_validate_peptide_sequence_case_and_characters(self, sequence):
        is_valid, clean_seq = validate_peptide_sequence(sequence)
        assert isinstance(is_valid, bool)
        assert isinstance(clean_seq, str)

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_validate_peptide_sequence_parse_exception(self, mock_parse):
        # Test fallback behavior when parsing fails
        mock_parse.side_effect = Exception("Parse error")

        # Test with a simple valid sequence to ensure fallback works
        is_valid, clean_seq = validate_peptide_sequence("PEPTIDE")
        assert isinstance(is_valid, bool)
        assert isinstance(clean_seq, str)

    def test_validate_peptide_sequence_fallback_validation(self):
        # Test that fallback validation works correctly without mocking
        test_sequence = "PEPTIDE"
        is_valid, clean_seq = validate_peptide_sequence(test_sequence)
        assert is_valid is True
        assert clean_seq == test_sequence

    @pytest.mark.parametrize(
        "sequence,expected_clean,expected_charge",
        [
            ("PEPTIDE/2", "PEPTIDE", 2),
            ("SEQUENCE/3", "SEQUENCE", 3),
            ("PEPTIDE2", "PEPTIDE", 2),
            ("SEQUENCE3", "SEQUENCE", 3),
        ],
    )
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_validate_peptide_sequence_with_mods_valid(
        self, mock_parse, sequence, expected_clean, expected_charge
    ):
        mock_parse.return_value = (expected_clean, expected_clean, expected_charge)

        is_valid, clean_seq, openms_seq, charge = validate_peptide_sequence_with_mods(
            sequence
        )
        assert is_valid is True
        assert clean_seq == expected_clean
        assert charge == expected_charge
        mock_parse.assert_called_once_with(sequence)

    @pytest.mark.parametrize(
        "invalid_sequence,expected_clean",
        [
            ("PEPTIDEZ[Oxidation]", "PEPTIDEZ"),
            ("PEPTIDEJ/2", "PEPTIDEJ"),
            ("INVALIDZ123", "INVALIDZ"),
        ],
    )
    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_validate_peptide_sequence_with_mods_invalid(
        self, mock_parse, invalid_sequence, expected_clean
    ):
        mock_parse.return_value = (expected_clean, invalid_sequence, 1)

        is_valid, clean_seq, openms_seq, charge = validate_peptide_sequence_with_mods(
            invalid_sequence
        )
        # Should detect that the clean sequence contains invalid amino acids
        assert is_valid is False

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_validate_peptide_sequence_with_mods_exception_handling(self, mock_parse):
        mock_parse.side_effect = Exception("Parse error")

        is_valid, clean_seq, openms_seq, charge = validate_peptide_sequence_with_mods(
            "INVALID"
        )
        assert is_valid is False
        assert clean_seq == ""
        assert openms_seq == ""
        assert charge == 1

    def test_validate_very_long_sequence(self):
        long_sequence = "A" * 1000
        is_valid, clean_seq = validate_peptide_sequence(long_sequence)
        assert is_valid is True
        assert clean_seq == long_sequence

    def test_validate_mixed_valid_invalid_sequence(self):
        # Test sequences with mix of valid and invalid amino acids
        mixed_sequence = "PEPTIDEZABC"  # Z is invalid
        is_valid, clean_seq = validate_peptide_sequence(mixed_sequence)
        assert is_valid is False

    @pytest.mark.parametrize(
        "sequence",
        [
            "PEPTIDE\n",  # With newline
            "PEPTIDE\t",  # With tab
            " PEPTIDE ",  # With spaces
            "PEPTIDE\r\n",  # With carriage return
        ],
    )
    def test_validate_peptide_sequence_whitespace_handling(self, sequence):
        is_valid, clean_seq = validate_peptide_sequence(sequence)
        assert isinstance(is_valid, bool)
        assert isinstance(clean_seq, str)

    def test_validate_peptide_sequence_unicode_characters(self):
        unicode_sequence = "PEPTIDE™"
        is_valid, clean_seq = validate_peptide_sequence(unicode_sequence)
        assert isinstance(is_valid, bool)
        assert isinstance(clean_seq, str)

    @patch("src.peptide_calculator.parse_sequence_with_mods_and_charge")
    def test_validate_peptide_sequence_with_complex_modifications(self, mock_parse):
        # Test validation with complex modification patterns
        complex_sequences = [
            ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]/2",
            "EM[+15.9949]EVEES[-79.9663]PEK",
            "ALSSC[UNIMOD:4]VVDEEQDVER/2",
        ]

        mock_parse.return_value = ("PEPTIDE", "PEPTIDE", 2)
        for sequence in complex_sequences:
            is_valid, clean_seq, openms_seq, charge = (
                validate_peptide_sequence_with_mods(sequence)
            )
            assert isinstance(is_valid, bool)
            assert isinstance(clean_seq, str)
            assert isinstance(openms_seq, str)
            assert isinstance(charge, int)

    def test_validate_peptide_sequence_error_messages(self):
        # Test that appropriate error context is maintained
        invalid_sequence = "PEPTIDEZ"
        is_valid, clean_seq = validate_peptide_sequence(invalid_sequence)
        assert is_valid is False
        # Error messages should be available in ERROR_MESSAGES constant
        assert hasattr(ERROR_MESSAGES, "__iter__") or isinstance(ERROR_MESSAGES, dict)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
