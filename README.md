<div align="center">
  <img src="assets/openms_transparent_bg_logo.svg" alt="OpenMS Logo" width="200" height="200" />
  <h3 align="center">OpenMS Peptide m/z Calculator</h3>
  <p align="center">An open-source tool for calculating peptide m/z values.</p>
</div>

<p align="center">
  <a href="https://github.com/OpenMS/mz-calculator/issues">Report Bug/Feature</a> •
  <a href="https://github.com/OpenMS/mz-calculator/pulls">Submit Code Changes</a>
</p>

# 📋 Overview

This web application provides accurate peptide mass-to-charge (m/z) ratio calculations for proteomics research.

## ✨ Key Features

- **Peptide m/z Calculations**: Calculate mass-to-charge ratios with support for modifications and charge states
- **Advanced Notation Support**:
  - Square bracket modifications: `M[Oxidation]PEPTIDE`
  - UNIMOD notation: `C[UNIMOD:4]PEPTIDE`
  - ProForma mass shifts: `PEPTIDE[+42.0106]`
  - Charge notation: `PEPTIDE/2`
- **Auto-Detection**: Automatically recognizes modifications and charge states from sequence notation

## 🚀 Quick Start

### Web Application

Access the calculator through your web browser - simply enter a peptide sequence and get instant m/z calculations.

### Example Usage

```
PEPTIDE                    # Basic sequence
M[Oxidation]PEPTIDE       # With modification
C[UNIMOD:4]PEPTIDE        # UNIMOD notation
PEPTIDE/2                 # With charge state
PEPTIDE[+42.0106]         # ProForma mass shift
```

## 💻 Installation

### Method: Direct Python Installation

This is the recommended method for development environments and direct Python execution.

#### Step 1: Create Python Virtual Environment

```bash
# Create virtual environment
python -m venv peptide-calculator-env

# Activate virtual environment
# On Windows:
peptide-calculator-env\Scripts\activate
# On macOS/Linux:
source peptide-calculator-env/bin/activate
```

#### Step 2: Install Dependencies

The application uses a compiled requirements file for reproducible installations:

```bash
# Install all dependencies
pip install -r requirements.txt

# Alternative: Install from pyproject.toml for development
pip install -e .
```

#### Step 3: Verify Core Dependencies

```bash
# Test PyOpenMS installation
python -c "import pyopenms; print('PyOpenMS version:', pyopenms.__version__)"

# Test Streamlit installation
streamlit --version
```

## 📖 Citation

If you use this work, please cite:

```
Müller, T. D., Siraj, A., et al. OpenMS WebApps: Building User-Friendly Solutions for MS Analysis. 
Journal of Proteome Research (2025). 
https://doi.org/10.1021/acs.jproteome.4c00872
```

## 📚 References

This work builds upon the following publications:

1. **Pfeuffer, J., Bielow, C., Wein, S. et al.** (2024). OpenMS 3 enables reproducible analysis of large-scale mass spectrometry data. *Nature Methods*, 21, 365–367. <https://doi.org/10.1038/s41592-024-02197-7>

2. **Röst, H. L., Schmitt, U., Aebersold, R., & Malmström, L.** (2014). pyOpenMS: a Python-based interface to the OpenMS mass-spectrometry algorithm library. *Proteomics*, 14(1), 74-77. <https://doi.org/10.1002/pmic.201300246> | PubMed: [24420968](https://pubmed.ncbi.nlm.nih.gov/24420968/)
