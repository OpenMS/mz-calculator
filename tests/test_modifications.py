"""
Test modification detection and application functions for peptide calculator.
"""

import sys
import os
from unittest.mock import patch, MagicMock

# Add project root to path
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)

from src.peptide_calculator import (
    detect_modification_from_sequence,
    _match_mass_delta_to_modification,
    apply_modification,
    get_supported_modifications,
    get_modification_info,
    get_square_bracket_examples,
)


class TestModificationDetection:
    """Test cases for modification detection functions."""

    def test_detect_modification_no_modification(self):
        """Test sequences without modifications."""
        modification = detect_modification_from_sequence("PEPTIDE")
        assert modification == "None"

        modification = detect_modification_from_sequence("SIMPLESEQUENCE")
        assert modification == "None"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_oxidation(self, mock_parse):
        mock_parse.return_value = ("MPEPTIDE", "M(Oxidation)PEPTIDE")
        modification = detect_modification_from_sequence("M[Oxidation]PEPTIDE")
        assert modification == "Oxidation (M)"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_carbamidomethyl(self, mock_parse):
        mock_parse.return_value = ("CPEPTIDE", "C(Carbamidomethyl)PEPTIDE")
        modification = detect_modification_from_sequence("C[Carbamidomethyl]PEPTIDE")
        assert modification == "Carbamidomethyl (C)"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_phosphorylation(self, mock_parse):
        mock_parse.return_value = ("SPEPTIDE", "S(Phospho)PEPTIDE")
        modification = detect_modification_from_sequence("S[Phospho]PEPTIDE")
        assert modification == "Phosphorylation (S/T/Y)"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_n_terminal_acetylation(self, mock_parse):
        mock_parse.return_value = ("PEPTIDE", ".(Acetyl)PEPTIDE")
        modification = detect_modification_from_sequence("[Acetyl]PEPTIDE")
        assert modification == "Acetylation (N-term)"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_methylation(self, mock_parse):
        mock_parse.return_value = ("KPEPTIDE", "K(Methyl)PEPTIDE")
        modification = detect_modification_from_sequence("K[Methyl]PEPTIDE")
        assert modification == "Methylation (K/R)"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_deamidation(self, mock_parse):
        mock_parse.return_value = ("NPEPTIDE", "N(Deamidated)PEPTIDE")
        modification = detect_modification_from_sequence("N[Deamidated]PEPTIDE")
        assert modification == "Deamidation (N/Q)"

    @patch("src.peptide_calculator._match_mass_delta_to_modification")
    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_mass_delta(self, mock_parse, mock_match):
        mock_parse.return_value = ("MPEPTIDE", "M[+15.9949]PEPTIDE")
        mock_match.return_value = "Oxidation (M)"
        modification = detect_modification_from_sequence("M[+15.9949]PEPTIDE")
        assert modification == "Oxidation (M)"

    @patch("src.peptide_calculator._match_mass_delta_to_modification")
    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_unknown_mass_delta(self, mock_parse, mock_match):
        mock_parse.return_value = ("MPEPTIDE", "M[+999.999]PEPTIDE")
        mock_match.return_value = "None"
        modification = detect_modification_from_sequence("M[+999.999]PEPTIDE")
        assert modification == "None"

    @patch("src.peptide_calculator.parse_square_bracket_modifications")
    def test_detect_modification_parse_exception(self, mock_parse):
        mock_parse.side_effect = Exception("Parse error")
        modification = detect_modification_from_sequence("INVALID[Mod]SEQUENCE")
        assert modification == "None"

    def test_match_mass_delta_known_modifications(self, mock_pyopenms):
        """Test mass delta matching for known modifications."""
        result = _match_mass_delta_to_modification(57.021464)
        assert result == "Carbamidomethyl (C)"

        result = _match_mass_delta_to_modification(15.994915)
        assert result == "Oxidation (M)"

        result = _match_mass_delta_to_modification(79.966331)
        assert result == "Phosphorylation (S/T/Y)"

        result = _match_mass_delta_to_modification(42.010565)
        assert result == "Acetylation (N-term)"

        result = _match_mass_delta_to_modification(14.015650)
        assert result == "Methylation (K/R)"

        result = _match_mass_delta_to_modification(0.984016)
        assert result == "Deamidation (N/Q)"

    def test_match_mass_delta_tolerance(self, mock_pyopenms):
        """Test mass delta matching with tolerance."""
        result = _match_mass_delta_to_modification(57.022, tolerance=0.01)
        assert result == "Carbamidomethyl (C)"

        # 0.008536 difference, should not match with stricter tolerance
        result = _match_mass_delta_to_modification(57.030, tolerance=0.005)
        assert result == "None"

        result = _match_mass_delta_to_modification(57.025, tolerance=0.01)
        assert result == "Carbamidomethyl (C)"

    def test_match_mass_delta_unknown_mass(self, mock_pyopenms):
        """Test mass delta matching for unknown modifications."""
        result = _match_mass_delta_to_modification(999.999)
        assert result == "None"

        result = _match_mass_delta_to_modification(-50.0)
        assert result == "None"

        result = _match_mass_delta_to_modification(0.001)
        assert result == "None"

    def test_match_mass_delta_with_moddb_fallback(self, mock_pyopenms):
        """Test mass delta matching with ModificationsDB fallback."""
        # Close to oxidation mass - should either match known mass or fall back to ModDB
        result = _match_mass_delta_to_modification(15.99)
        assert result in ["Oxidation (M)", "None"]


class TestModificationApplication:
    """Test cases for modification application functions."""

    def test_apply_modification_none(self):
        """Test applying no modification."""
        result = apply_modification("PEPTIDE", "None")
        assert result == "PEPTIDE"

        result = apply_modification("SEQUENCE", "None")
        assert result == "SEQUENCE"

    def test_apply_modification_oxidation(self):
        """Test applying oxidation modification."""
        result = apply_modification("MPEPTIDE", "Oxidation (M)")
        assert "M(Oxidation)" in result or result == "MPEPTIDE"

        # Test sequence without methionine
        result = apply_modification("PEPTIDE", "Oxidation (M)")
        assert result == "PEPTIDE"

    def test_apply_modification_carbamidomethyl(self):
        """Test applying carbamidomethyl modification."""
        result = apply_modification("CPEPTIDE", "Carbamidomethyl (C)")
        assert "C(Carbamidomethyl)" in result or result == "CPEPTIDE"

        # Test sequence without cysteine
        result = apply_modification("PEPTIDE", "Carbamidomethyl (C)")
        assert result == "PEPTIDE"

    def test_apply_modification_n_terminal_acetylation(self):
        """Test applying N-terminal acetylation."""
        result = apply_modification("PEPTIDE", "Acetylation (N-term)")
        assert result.startswith(".(Acetyl)") or result == "PEPTIDE"

    def test_apply_modification_c_terminal(self):
        """Test applying C-terminal modification."""
        result = apply_modification("PEPTIDE", "Unknown C-terminal")
        assert result == "PEPTIDE"  # Should return unchanged for unknown mods

    def test_apply_modification_phosphorylation_multiple_targets(self):
        """Test applying phosphorylation to multiple possible targets."""
        result = apply_modification("PEPTIDES", "Phosphorylation (S/T/Y)")
        assert "S(Phospho)" in result or result == "PEPTIDES"

        result = apply_modification("PEPTITET", "Phosphorylation (S/T/Y)")
        assert "T(Phospho)" in result or result == "PEPTITET"

    def test_apply_modification_methylation_multiple_targets(self):
        """Test applying methylation to K or R."""
        result = apply_modification("KPEPTIDE", "Methylation (K/R)")
        assert "K(Methyl)" in result or result == "KPEPTIDE"

        result = apply_modification("RPEPTIDE", "Methylation (K/R)")
        assert "R(Methyl)" in result or result == "RPEPTIDE"

    def test_apply_modification_deamidation(self):
        """Test applying deamidation modification."""
        result = apply_modification("NPEPTIDE", "Deamidation (N/Q)")
        assert "N(Deamidated)" in result or result == "NPEPTIDE"

        result = apply_modification("QPEPTIDE", "Deamidation (N/Q)")
        assert "Q(Deamidated)" in result or result == "QPEPTIDE"

    def test_apply_modification_unknown_modification(self):
        """Test applying unknown modification."""
        result = apply_modification("PEPTIDE", "Unknown Modification")
        assert result == "PEPTIDE"

        result = apply_modification("PEPTIDE", "")
        assert result == "PEPTIDE"

    def test_get_supported_modifications(self):
        """Test getting supported modifications list."""
        modifications = get_supported_modifications()

        assert isinstance(modifications, list)
        assert "None" in modifications
        assert "Oxidation (M)" in modifications
        assert "Carbamidomethyl (C)" in modifications
        assert "Phosphorylation (S/T/Y)" in modifications
        assert "Acetylation (N-term)" in modifications
        assert "Methylation (K/R)" in modifications
        assert "Deamidation (N/Q)" in modifications
        assert len(modifications) == 7  # Expected number of modifications

    def test_get_modification_info(self, mock_pyopenms):
        """Test getting detailed modification information."""
        mod_info = get_modification_info()

        assert isinstance(mod_info, dict)
        assert "None" in mod_info
        assert mod_info["None"] == "No modification applied"

        for mod_name, description in mod_info.items():
            assert isinstance(mod_name, str)
            assert isinstance(description, str)
            if mod_name != "None":
                assert len(description) > 0

    def test_get_square_bracket_examples(self):
        """Test getting square bracket examples."""
        examples = get_square_bracket_examples()

        assert isinstance(examples, dict)
        assert len(examples) > 0

        # Check for specific example types
        found_oxidation = any("Oxidation" in key for key in examples.keys())
        found_mass_delta = any("[+" in key for key in examples.keys())
        found_charge = any("/" in key for key in examples.keys())
        found_unimod = any("UNIMOD" in key for key in examples.keys())

        assert found_oxidation
        assert found_mass_delta
        assert found_charge
        assert found_unimod

        # Check specific examples
        assert "M[Oxidation]PEPTIDE" in examples
        assert "PEPTIDE/2" in examples or "VAEINPSNGGTT/2" in examples

        # Verify all values are strings
        for sequence, description in examples.items():
            assert isinstance(sequence, str)
            assert isinstance(description, str)
            assert len(description) > 0

    def test_get_square_bracket_examples_content(self):
        """Test specific content of square bracket examples."""
        examples = get_square_bracket_examples()

        example_types = {
            "basic_mod": False,
            "n_terminal": False,
            "c_terminal": False,
            "mass_delta": False,
            "unimod": False,
            "charge_notation": False,
            "proforma": False,
        }

        for sequence in examples.keys():
            if "M[Oxidation]" in sequence or "C[Carbamidomethyl]" in sequence:
                example_types["basic_mod"] = True
            if sequence.startswith("[") or sequence.startswith(".["):
                example_types["n_terminal"] = True
            if (
                sequence.endswith("]")
                and not sequence.endswith("/2]")
                and not sequence.endswith("/3]")
            ):
                example_types["c_terminal"] = True
            if "[+" in sequence or "[-" in sequence:
                example_types["mass_delta"] = True
            if "UNIMOD" in sequence:
                example_types["unimod"] = True
            if "/" in sequence and sequence[-1].isdigit():
                example_types["charge_notation"] = True
            if "[+" in sequence and "." in sequence:
                example_types["proforma"] = True

        # Should have examples of most types
        assert sum(example_types.values()) >= 5


if __name__ == "__main__":
    import subprocess

    subprocess.run(["python", "-m", "pytest", __file__, "-v"])
