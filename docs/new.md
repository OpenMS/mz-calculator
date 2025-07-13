# New Features and Updates

## Latest Updates

### Version 1.0 - Enhanced Peptide Calculator

#### 🚀 New Features

**ProForma Arbitrary Mass Shifts Support**
- Support for arbitrary mass deltas in ProForma format
- Examples: `PEPTIDE[+15.9949]`, `SEQUENCE[-18.0106]`
- Scientific notation support: `PEPTIDE[+1.5e2]`
- Multiple mass shifts in single sequence: `EM[+15.9949]EVEES[-79.9663]PEK`

**Enhanced UNIMOD Integration**
- Full UNIMOD notation support: `C[UNIMOD:4]`, `M[UNIMOD:35]`
- Case-insensitive UNIMOD IDs: `[unimod:4]` = `[UNIMOD:4]`
- Automatic UNIMOD to mass delta conversion
- Standardized modification references

**Improved Auto-Detection**
- Real-time modification detection from sequence notation
- Automatic charge state extraction from sequences
- Smart dropdown updates based on detected modifications
- Synchronized input fields for better user experience

**Advanced Charge Notation**
- Slash format: `PEPTIDE/2`, `SEQUENCE/3`
- Trailing number format: `PEPTIDE2`, `SEQUENCE3`
- Combined with modifications: `M[Oxidation]PEPTIDE/2`
- Auto-update of charge input field

#### 🔧 Technical Improvements

**Enhanced Parser Engine**
- Direct pyOpenMS parsing for ProForma sequences
- Fallback parsing for complex notation
- Better error handling and validation
- Support for leading dot notation: `.PEPTIDE`

**Robust Modification Handling**
- Lazy loading of pyOpenMS ModificationsDB
- Cached modification lookups for performance
- Mass delta to modification name matching
- Support for N-terminal and C-terminal modifications

**Improved Validation**
- Comprehensive amino acid validation
- Real-time sequence analysis
- Better error messages with troubleshooting tips
- Invalid character detection and reporting

#### 🎯 User Experience Enhancements

**Smart Input Processing**
- Auto-completion of modification dropdown
- Real-time charge state detection
- Session state management for better performance
- Caching of expensive operations

**Enhanced Results Display**
- Detailed sequence information
- Charge source indication (🔗 for auto-detected)
- Amino acid composition breakdown
- Expandable additional information

**Better Documentation**
- Comprehensive examples in UI
- Interactive help sections
- Troubleshooting guides
- Advanced notation examples

#### 📊 Performance Optimizations

**Caching Improvements**
- Session-based analysis caching
- Modification list caching (1-hour TTL)
- Example data caching
- Reduced redundant calculations

**Memory Management**
- Lazy loading of heavy dependencies
- Cached ModificationsDB instance
- Optimized sequence parsing
- Reduced memory footprint

## Previous Updates

### Version 0.5 - Core Calculator

**Initial Features**
- Basic peptide m/z calculation
- Common modification support
- Charge state handling
- pyOpenMS integration

**Supported Modifications**
- Oxidation (M)
- Carbamidomethyl (C)
- Phosphorylation (S/T/Y)
- Acetylation (N-term)
- Methylation (K/R)
- Deamidation (N/Q)

### Version 0.1 - Foundation

**Core Functionality**
- Streamlit-based web interface
- Basic sequence input
- Simple m/z calculations
- Result display

