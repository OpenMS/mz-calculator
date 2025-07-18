# User Guide

## Overview

The Peptide M/Z Calculator is a web-based application that calculates accurate mass-to-charge (m/z) ratios for peptides used in mass spectrometry and proteomics research. Built with pyOpenMS, it provides precise calculations for peptides with various modifications and charge states.

## Getting Started

### Basic Usage

1. **Enter Peptide Sequence**: Type your peptide sequence using standard amino acid codes (e.g., `PEPTIDE`)
2. **Select Modifications** (Optional): Choose from common modifications or use notation
3. **Set Charge State**: Specify the charge state (default is 2)
4. **Calculate**: Click "Calculate m/z" to get results

### Supported Amino Acids

The calculator supports all standard amino acids:
`A, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V, X, U`

- **X**: Any amino acid (unknown)
- **U**: Selenocysteine

## Input Formats

### 1. Basic Sequences
```
PEPTIDE
STREAMLIT
ALSSCVVDEEQDVER
```

### 2. Modification Notation

#### Square Bracket Modifications
```
M[Oxidation]PEPTIDE          # Methionine oxidation
C[Carbamidomethyl]PEPTIDE    # Carbamidomethylated cysteine
[Acetyl]PEPTIDE              # N-terminal acetylation
PEPTIDE[Amidated]            # C-terminal amidation
```

#### UNIMOD Notation (Standardized)
```
M[UNIMOD:35]PEPTIDE          # UNIMOD:35 = Oxidation
C[UNIMOD:4]PEPTIDE           # UNIMOD:4 = Carbamidomethyl
[UNIMOD:1]PEPTIDE            # UNIMOD:1 = Acetyl (N-term)
```

#### ProForma Arbitrary Mass Shifts
```
PEPTIDE[+15.9949]            # +15.9949 Da mass shift
EM[+15.9949]EVEES[-79.9663]PEK  # Multiple mass shifts
PEPTIDE[+1.5e2]              # Scientific notation
SEQUENCE[-18.0106]           # Negative mass shift
```

### 3. Charge Notation

#### Slash Format
```
PEPTIDE/2                    # Charge state 2
M[Oxidation]PEPTIDE/3        # Modified peptide, charge 3
```

#### Trailing Number Format
```
PEPTIDE2                     # Charge state 2
SEQUENCE3                    # Charge state 3
```

### 4. Combined Notation
```
.[Acetyl]M[Oxidation]PEPTIDE/2      # N-term acetyl + Met oxidation + charge 2
QVVPC[+57.021464]STSER/2            # Mass delta + charge notation
ALSSC[UNIMOD:4]VVDEEQDVER/2         # UNIMOD + charge notation
```

## Modifications

### Common Modifications Available

| Modification | Target | Mass Delta (Da) |
|--------------|--------|----------------|
| Oxidation | M | +15.994915 |
| Carbamidomethyl | C | +57.021464 |
| Phosphorylation | S/T/Y | +79.966331 |
| Acetylation | N-term | +42.010565 |
| Methylation | K/R | +14.015650 |
| Deamidation | N/Q | +0.984016 |

### Auto-Detection Features

- **Modification Detection**: The calculator automatically detects modifications from sequence notation
- **Charge Detection**: Charge states are automatically extracted from notation
- **Dropdown Updates**: The modification dropdown updates when modifications are detected
- **Input Synchronization**: Charge input field updates when charge notation is found

## Results Interpretation

### Primary Results
- **m/z Ratio**: The calculated mass-to-charge ratio
- **Monoisotopic Mass**: The exact mass in Daltons
- **Charge State**: Final charge state used (with source indicator)

### Sequence Information
- **Stripped Sequence**: Clean amino acid sequence
- **Modified Sequence**: Sequence with applied modifications
- **Molecular Formula**: Chemical formula of the peptide

### Additional Information
- **Sequence Length**: Number of amino acids
- **Applied Modification**: Details of modifications used
- **Charge Source**: Whether charge came from notation or input
- **Amino Acid Composition**: Count of each amino acid

## Advanced Features

### Leading Dot Notation
```
.LLVLPKFGM[+15.9949]LMLGPDDFR
```
The leading dot (.) represents the N-terminus and is automatically handled.

### Multiple Modifications
```
.[Acetyl]M[Oxidation]PEPTIDEC[Carbamidomethyl]
```
Multiple modifications can be applied to different positions.

### Terminal Modifications
```
[Acetyl]PEPTIDE              # N-terminal
PEPTIDE.[Amidated]           # C-terminal (dot required)
```

## Troubleshooting

### Common Issues

1. **Invalid Amino Acid Error**
   - Check for typos in sequence
   - Ensure only valid amino acid codes are used
   - Remove any special characters

2. **Modification Not Recognized**
   - Try using UNIMOD notation
   - Use mass delta format: `[+mass]`
   - Check modification spelling

3. **Charge State Issues**
   - Ensure charge is between 1-20
   - Check charge notation format
   - Verify slash or trailing number format

4. **Sequence Parsing Errors**
   - Try without modifications first
   - Check bracket pairing `[]`
   - Verify modification placement

### Best Practices

1. **Start Simple**: Test with basic sequences first
2. **Use Standard Notation**: Prefer UNIMOD IDs for consistency
3. **Validate Input**: Check sequences before complex modifications
4. **Check Results**: Verify m/z values against expected ranges
5. **Use Auto-Detection**: Let the calculator detect modifications and charges

## Tips for Efficient Use

- Use the auto-detection feature by entering complete notation
- The modification dropdown updates automatically
- Charge fields synchronize with sequence notation
- ProForma arbitrary mass shifts provide maximum flexibility
- UNIMOD notation ensures standardized modifications
- Combine different notation types as needed

## Technical Details

### Calculation Method
- Uses pyOpenMS AASequence class for accuracy
- Monoisotopic masses for all calculations
- m/z = (Monoisotopic Mass + Charge × Proton Mass) / Charge
- Proton mass: 1.007276466879 Da

### Supported Formats
- ProForma notation (partial support)
- UNIMOD standard modifications
- OpenMS native format
- Custom mass delta notation
- Multiple charge notation formats