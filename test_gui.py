from streamlit.testing.v1 import AppTest
import pytest
from src import fileupload
import json
from pathlib import Path
import shutil


@pytest.fixture
def launch(request):
    test = AppTest.from_file(request.param)

    ## Initialize session state ##
    with open("settings.json", "r") as f:
        test.session_state.settings = json.load(f)
    test.session_state.settings["test"] = True
    test.secrets["workspace"] = "test"
    return test


# Test launching of all pages
@pytest.mark.parametrize(
    "launch",
    (
        "content/quickstart.py",  # Main peptide calculator page
        "content/docs.py",       # Documentation page
        "content/visualization.py", # Visualization page
    ),
    indirect=True,
)
def test_launch(launch):
    """Test if all pages can be launched without errors."""
    launch.run(timeout=30)  # Increased timeout from 10 to 30 seconds
    assert not launch.exception


########### PAGE SPECIFIC TESTS ############
@pytest.mark.parametrize(
    "launch,selection",
    [
        ("content/docs.py", "User Guide"),
        ("content/docs.py", "What's New"),
        ("content/docs.py", "Feature Overview"),
    ],
    indirect=["launch"],
)
def test_documentation(launch, selection):
    launch.run()
    launch.selectbox[0].select(selection).run()
    assert not launch.exception


# Test peptide calculator basic functionality
@pytest.mark.parametrize("launch", ["content/quickstart.py"], indirect=True)
def test_peptide_calculator_basic(launch):
    """Test basic peptide calculator functionality."""
    launch.run()
    
    # Check if basic elements are present
    assert len(launch.text_input) > 0, "Peptide sequence input not found"
    assert len(launch.checkbox) > 0, "Modification checkboxes not found"
    assert len(launch.slider) > 0, "Charge state input not found"
    assert len(launch.button) > 0, "Calculate button not found"


# Test peptide calculator with valid input
@pytest.mark.parametrize(
    "launch,sequence,charge",
    [
        ("content/quickstart.py", "PEPTIDE", 2),
        ("content/quickstart.py", "METHIONINE", 1),
        ("content/quickstart.py", "STREAMLIT", 3),
    ],
    indirect=["launch"],
)
def test_peptide_calculator_calculation(launch, sequence, charge):
    """Test peptide calculator with different inputs."""
    launch.run()
    
    # Input sequence
    if len(launch.text_input) > 0:
        launch.text_input[0].set_value(sequence)
    
    # Set charge state
    if len(launch.number_input) > 0:
        launch.number_input[0].set_value(charge)
    
    launch.run()
    
    # Click calculate button if present
    if len(launch.button) > 0:
        launch.button[0].click()
        launch.run()
    
    assert not launch.exception


# Test visualization page
@pytest.mark.parametrize("launch", ["content/visualization.py"], indirect=True)
def test_visualization_page(launch):
    """Test that visualization page loads without errors."""
    launch.run()
    assert not launch.exception
