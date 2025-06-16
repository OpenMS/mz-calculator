import pytest
import sys
import os
from unittest.mock import MagicMock, patch

# Add project root to path for imports
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, PROJECT_ROOT)


def pytest_configure(config):
    """
    Configure pytest with custom markers for test categorization.
    These markers help categorize and run specific sets of tests.
    """
    config.addinivalue_line("markers", "unit: Unit tests for individual functions")
    config.addinivalue_line("markers", "integration: Integration tests")
    config.addinivalue_line("markers", "validation: Validation function tests")
    config.addinivalue_line("markers", "parsing: Parsing function tests")
    config.addinivalue_line("markers", "modification: Modification handling tests")
    config.addinivalue_line("markers", "calculation: M/Z calculation tests")
    config.addinivalue_line("markers", "edge_case: Edge case tests")


@pytest.fixture(scope="session")
def mock_pyopenms():
    """
    Mocks the pyOpenMS library to enable testing of
    peptide m/z calculation logic without requiring
    the actual pyOpenMS dependency to be installed.
    This fixture provides mock objects for AASequence and ModificationsDB
    with predefined return values for common methods, simulating
    the behavior of pyOpenMS for testing purposes.
    """
    mock_poms = MagicMock()
    mock_poms.__version__ = "2.9.1"

    mock_aa_seq = MagicMock()
    mock_aa_seq.getMZ.return_value = 500.25
    mock_aa_seq.getMonoWeight.return_value = 999.49
    mock_aa_seq.getFormula.return_value.toString.return_value = "C43H66N12O12"
    mock_poms.AASequence.fromString.return_value = mock_aa_seq

    mock_mod_db = MagicMock()
    mock_mod = MagicMock()
    mock_mod.getId.return_value = "Oxidation"
    mock_mod.getFullName.return_value = "Oxidation"
    mock_mod.getDiffMonoMass.return_value = 15.994915
    mock_mod.getTermSpecificity.return_value = 0  # ANYWHERE
    mock_mod.getOrigin.return_value = "M"
    mock_mod_db.getModification.return_value = mock_mod
    mock_poms.ModificationsDB.return_value = mock_mod_db

    mock_poms.ResidueModification.ANYWHERE = 0
    mock_poms.ResidueModification.N_TERM = 1
    mock_poms.ResidueModification.C_TERM = 2

    with patch("pyopenms", mock_poms):
        with patch.dict(sys.modules, {"pyopenms": mock_poms}):
            yield mock_poms


@pytest.fixture
def sample_peptides():
    """
    Provides a dictionary of sample peptide sequences categorized
    by their characteristics (e.g., basic, with modifications, with charge, invalid).
    Useful for testing various peptide parsing and calculation scenarios.
    """
    return {
        "basic": ["PEPTIDE", "SEQUENCE", "PROTEIN", "ANALYSIS"],
        "with_modifications": [
            "M[Oxidation]PEPTIDE",
            "C[Carbamidomethyl]PEPTIDE",
            "[Acetyl]PEPTIDE",
            "PEPTIDE[Amidated]",
        ],
        "with_charge": ["PEPTIDE/2", "SEQUENCE/3", "PEPTIDE2", "SEQUENCE3"],
        "mass_deltas": [
            "M[+15.9949]PEPTIDE",
            "C[+57.021464]PEPTIDE",
            "PEPTIDE[+42.0106]",
        ],
        "complex": [
            ".[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]/2",
            "EM[+15.9949]EVEES[-79.9663]PEK",
            "ALSSC[UNIMOD:4]VVDEEQDVER/2",
        ],
        "invalid": ["PEPTIDEZ", "PEPTIDE123", "", "   "],
    }


@pytest.fixture
def expected_modifications():
    """
    Provides a dictionary of common peptide modifications and their
    expected monoisotopic mass deltas.
    Useful for validating modification parsing and mass calculations.
    """
    return {
        "Oxidation (M)": 15.994915,
        "Carbamidomethyl (C)": 57.021464,
        "Phosphorylation (S/T/Y)": 79.966331,
        "Acetylation (N-term)": 42.010565,
        "Methylation (K/R)": 14.015650,
        "Deamidation (N/Q)": 0.984016,
    }
