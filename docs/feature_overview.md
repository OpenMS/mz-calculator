# Features Overview

## Core Functionality

### 🧬 Peptide Mass Calculation Engine

**Accurate m/z Calculations**
- Powered by pyOpenMS for precision mass spectrometry calculations
- Monoisotopic mass calculations using exact atomic masses
- Support for charge states from 1 to 20
- Professional-grade accuracy for research applications

**Advanced Sequence Processing**
- Real-time sequence validation and analysis
- Support for all standard amino acids (20 + X, U)
- Intelligent error detection and reporting
- Comprehensive sequence length validation

### 🔬 Modification Support

**Multiple Notation Formats**
- **Square Bracket Notation**: `M[Oxidation]PEPTIDE`
- **UNIMOD Standard**: `C[UNIMOD:4]PEPTIDE`
- **ProForma Mass Deltas**: `PEPTIDE[+15.9949]`
- **OpenMS Native Format**: `PEPTIDE(Oxidation)`

**Comprehensive Modification Database**
- 6+ common modifications pre-configured
- Full UNIMOD database integration
- Custom mass delta support
- N-terminal and C-terminal modifications

**Supported Common Modifications**
- Oxidation (M): +15.994915 Da
- Carbamidomethyl (C): +57.021464 Da
- Phosphorylation (S/T/Y): +79.966331 Da
- Acetylation (N-term): +42.010565 Da
- Methylation (K/R): +14.015650 Da
- Deamidation (N/Q): +0.984016 Da

### ⚡ Smart Auto-Detection

**Real-Time Analysis**
- Automatic modification detection from sequence notation
- Charge state extraction from sequence format
- Smart dropdown updates based on detected features
- Synchronized input fields for seamless user experience

**Intelligent Parsing**
- Direct pyOpenMS parsing for ProForma sequences
- Fallback parsing for complex notation combinations
- Leading dot notation handling (N-terminus indicator)
- Multiple modification support in single sequence

### 📊 Performance Optimizations

**Caching System**
- Session-based sequence analysis caching
- Modification list caching (1-hour TTL)
- Example data caching for improved load times
- Reduced redundant pyOpenMS operations

**Memory Management**
- Lazy loading of pyOpenMS modules
- Cached ModificationsDB instance
- Optimized sequence parsing algorithms
- Minimal memory footprint design

## User Interface Features

### 🎨 Modern Web Interface

**Responsive Design**
- Clean, professional layout
- Mobile-friendly responsive design
- Intuitive navigation and controls
- Streamlit-powered interactive components

**Smart Input Processing**
- Auto-complete sequence suggestions
- Real-time validation feedback
- Context-sensitive help and examples
- Progress indicators for calculations

### 📈 Comprehensive Results Display

**Primary Results**
- **m/z Ratio**: Precise mass-to-charge calculation
- **Monoisotopic Mass**: Exact mass in Daltons
- **Molecular Formula**: Complete chemical formula
- **Charge State**: With source indication (🔗 for auto-detected)

**Detailed Analysis**
- Original vs. modified sequence comparison
- Amino acid composition breakdown
- Sequence length and statistics
- Applied modification details

**Additional Information**
- Expandable sections for advanced data
- Modification source tracking
- Charge state origin indication
- Comprehensive troubleshooting guidance

### 📚 Interactive Documentation

**In-App Help**
- Comprehensive notation examples
- Interactive expandable sections
- Real-time help contextual to input
- Advanced notation tutorials

**Example Gallery**
- 20+ example sequences with descriptions
- ProForma notation examples
- UNIMOD reference examples
- Charge notation demonstrations

## Advanced Features

### 🧪 ProForma Support

**Arbitrary Mass Shifts**
- Support for any mass delta: `[+42.0106]`, `[-18.0106]`
- Scientific notation: `[+1.5e2]`, `[-2.5e-1]`
- Multiple mass shifts per sequence
- Flexible positioning (any amino acid)

**ProForma Compliance**
- PSI-MS ProForma notation standard
- Direct pyOpenMS parsing when possible
- Fallback conversion for complex cases
- Future-proof design for standard evolution

### 🔭 UNIMOD Integration

**Standardized Modifications**
- Full UNIMOD ID support: `[UNIMOD:4]`, `[UNIMOD:35]`
- Case-insensitive notation: `[unimod:4]` = `[UNIMOD:4]`
- Automatic mass delta resolution
- Standardized modification names

**Database Features**
- Cached UNIMOD database access
- Fast modification lookup
- Mass-to-modification matching
- Tolerance-based identification

### 🔄 Multi-Format Charge Notation

**Flexible Input Formats**
- **Slash notation**: `PEPTIDE/2`, `SEQUENCE/3`
- **Trailing numbers**: `PEPTIDE2`, `SEQUENCE3`
- **Combined with mods**: `M[Oxidation]PEPTIDE/2`
- **Auto-detection**: Input field synchronization

**Intelligent Processing**
- Charge state validation (1-20 range)
- Format recognition and parsing
- Override protection for user input
- Clear source indication in results

## Technical Architecture

### 🔧 Backend Engine

**pyOpenMS Integration**
- AASequence class for calculations
- ModificationsDB for modification data
- Formula class for molecular formulas
- Exact mass calculations using isotope data

**Error Handling**
- Comprehensive exception catching
- User-friendly error messages
- Detailed troubleshooting guidance
- Graceful fallback mechanisms

### 📊 Data Processing Pipeline

**Sequence Analysis Flow**
1. Input validation and sanitization
2. Charge notation extraction
3. Modification parsing and conversion
4. pyOpenMS sequence object creation
5. Mass and formula calculation
6. Results formatting and display

**Caching Strategy**
- Function-level caching for expensive operations
- Session state management
- TTL-based cache expiration
- Memory-conscious cache sizing

## Quality Assurance

### 🧪 Testing Framework

**Comprehensive Test Suite**
- Unit tests for all calculation functions
- Integration tests for complete workflows
- Validation tests for input parsing
- Modification handling tests
- Edge case and error condition tests

**Continuous Integration**
- Automated testing on code changes
- Multiple Python version compatibility
- Cross-platform testing
- Performance regression testing

### 🔒 Data Accuracy

**Validation Systems**
- Amino acid sequence validation
- Modification syntax checking
- Charge state range validation
- Mass calculation verification

**Reference Standards**
- UNIMOD standard modification database
- PSI-MS ProForma notation compliance
- OpenMS calculation accuracy
- Mass spectrometry community standards

## Performance Metrics

### 📈 Scalability Features

- Session-based state management
- Efficient caching strategies
- Lazy loading for heavy dependencies
- Optimized rendering pipeline