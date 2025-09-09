# Gene Annotation Analysis

A scientific gene annotation application that analyzes BLAST results and provides comprehensive gene information using the Ensembl REST API.

## 🧬 Project Overview

The Gene Annotation Analysis tool is a rule-compliant scientific application that:
- **Processes BLAST Results**: Reads CSV files from BLAST searches
- **Gene Annotation**: Integrates with Ensembl API for gene information
- **Intelligent Column Detection**: Automatically finds required columns
- **Data Validation**: Ensures scientific data integrity
- **Structured Error Handling**: Provides actionable error responses with citations

## 🏗️ Architecture

```
Genome_Analyzer_Data/
├── src/                    # Source code
│   ├── core/              # Core analysis engine
│   │   ├── gene_annotator.py    # Ensembl API integration
│   │   └── data_processor.py    # CSV processing and validation
│   ├── utils/             # Utility functions
│   │   └── error_handler.py     # Structured error handling
│   └── config/            # Configuration files
│       └── constants.py         # Scientific constants and thresholds
├── gene_annotation_script.py    # Main application script
├── requirements.txt              # Python dependencies
└── README.md                    # This file
```

## 🚀 Features

### Core Capabilities
- **CSV File Processing**: Load and validate BLAST result files
- **Column Detection**: Intelligent detection of required columns
- **Gene Annotation**: Ensembl API integration for gene information
- **Data Validation**: Comprehensive coordinate and data validation
- **Progress Tracking**: Real-time progress updates during processing

### Scientific Validation
- **Real Data Only**: No mock data or fallbacks
- **API Compliance**: Proper Ensembl API integration with rate limiting
- **Error Traceability**: Complete error context and scientific citations
- **Coordinate Validation**: Genomic coordinate integrity checks

### Performance Features
- **Rate Limiting**: Respectful API usage (15 requests/second)
- **Progress Reporting**: Updates every 10 rows
- **Batch Processing**: Efficient handling of large datasets
- **Memory Management**: Automatic resource cleanup

## 📋 Requirements

### System Requirements
- Python 3.8+
- Windows 10/11 (PowerShell compatible)
- Internet connection for Ensembl API access
- 512MB+ RAM for large datasets

### Python Dependencies
```bash
pip install -r requirements.txt
```

**Core Dependencies:**
- `requests>=2.31.0` - HTTP API communication
- `pandas>=2.0.0` - Data manipulation and CSV processing
- `numpy>=1.24.0` - Numerical computing

## 🛠️ Installation

1. **Download** the project files
2. **Install Dependencies**:
   ```bash
   cd Genome_Analyzer_Data
   pip install -r requirements.txt
   ```
3. **Verify Installation**:
   ```bash
   python gene_annotation_script.py --test-api
   ```

## 📖 Usage

### Basic Usage

Process a single CSV file:
```bash
python gene_annotation_script.py --file blast_results.csv
```

### Advanced Usage

Process multiple files with custom output:
```bash
python gene_annotation_script.py --file file1.csv file2.csv --output results/
```

Process all CSV files in a directory:
```bash
python gene_annotation_script.py --directory data/ --pattern "*.csv"
```

### Command Line Options

| Option | Description | Required |
|--------|-------------|----------|
| `--file` | CSV file(s) to process | Yes* |
| `--directory` | Directory containing CSV files | Yes* |
| `--pattern` | File pattern for directory processing | No |
| `--output` | Output directory for annotated files | No |
| `--columns` | Required column types | No |
| `--log-level` | Logging level (DEBUG/INFO/WARNING/ERROR) | No |
| `--test-api` | Test Ensembl API connection | No |

*Either `--file` or `--directory` is required

## 🔬 Scientific Parameters

### Gene Annotation Parameters
- **API Server**: Ensembl REST API (https://rest.ensembl.org)
- **Rate Limit**: 15 requests per second
- **Timeout**: 30 seconds for gene annotation requests
- **Coordinate Validation**: 0 to 250,000,000 bp range

### Column Detection
The tool automatically detects common column name variations:

**Accession Columns:**
- `subject acc.ver`, `accession`, `acc`, `subject_accession`

**Start Coordinate Columns:**
- `s. start`, `start`, `start_pos`, `start_position`

**End Coordinate Columns:**
- `s. end`, `s. End`, `Subject End`, `end`, `end_pos`

## 📊 Output Format

### Annotated CSV Files
The tool creates new CSV files with the `_annotated` suffix containing:

- All original columns from the input file
- **Gene_Annotation** column with gene information:
  - `[GeneName](in)` - Sequence completely within gene
  - `[GeneName]Overlap` - Sequence partially overlaps gene
  - `[Gene1]-------[Gene2]` - Sequence overlaps multiple genes
  - `Intergenic` - No gene overlap
  - `Non-standard Chromosome/Contig` - Invalid chromosome format
  - `Invalid Coordinates` - Coordinate validation failure

### Processing Summary
After completion, the tool displays:
- Total files processed
- Success/failure counts
- Total rows processed
- Output file locations

## 🧪 Testing

### Test API Connection
```bash
python gene_annotation_script.py --test-api
```

### Run with Debug Logging
```bash
python gene_annotation_script.py --file test.csv --log-level DEBUG
```

## 📝 Logging

### Log Files
- `logs/gene_annotation.log` - General application logs
- `logs/error.log` - Error-specific logs

### Log Levels
- **INFO**: General operations and progress
- **WARNING**: Non-critical issues
- **ERROR**: Errors with structured responses
- **DEBUG**: Detailed debugging information

## 🔧 Configuration

### Scientific Constants
All scientific parameters are defined in `src/config/constants.py` with:
- Scientific citations and references
- Validated threshold values
- API rate limiting parameters
- Coordinate validation settings

### Customization
Modify constants in `src/config/constants.py` for:
- Different API endpoints
- Custom rate limiting
- Coordinate range adjustments
- Column name variations

## 🚨 Error Handling

### Structured Error Responses
All errors return structured information with:
- **type**: Error category
- **message**: Human-readable description
- **citation**: Scientific reference
- **recommended_next_step**: Actionable guidance

### Error Types
- **API_ERROR**: External API communication issues
- **DATA_VALIDATION_ERROR**: Input validation failures
- **GENE_ANNOTATION_ERROR**: Gene annotation failures
- **FILE_PROCESSING_ERROR**: File processing issues

## 🔄 Data Processing Pipeline

### 1. File Discovery
- Find CSV files based on input parameters
- Validate file existence and accessibility

### 2. Column Detection
- Intelligently detect required columns
- Handle common column name variations
- Validate column presence and types

### 3. Data Validation
- Convert coordinates to numeric format
- Validate coordinate ranges and relationships
- Clean invalid or missing data

### 4. Gene Annotation
- Extract chromosome information from accessions
- Query Ensembl API for gene information
- Process and format gene overlap data

### 5. Output Generation
- Create annotated CSV files
- Generate processing summaries
- Log all operations and results

## 🔬 Scientific Validation

### Data Sources
- **Ensembl API**: Gene annotation data (GRCh38 assembly)
- **NCBI Accessions**: Chromosome identification
- **Coordinate System**: Ensembl coordinate system

### Quality Assurance
- **No Mock Data**: All data must be scientifically validated
- **API Compliance**: Proper external service integration
- **Error Traceability**: Complete error context and citations
- **Coordinate Validation**: Genomic coordinate integrity

## 🚀 Future Enhancements

### Planned Features
- **Batch Processing**: Process multiple files simultaneously
- **Caching System**: Cache gene annotations for performance
- **Advanced Visualization**: Interactive gene browser
- **Multi-Organism Support**: Beyond human genome analysis

### API Expansions
- **Additional Databases**: Expanded gene annotation sources
- **Custom Algorithms**: User-defined annotation methods
- **Performance Optimization**: Parallel processing capabilities

## 📚 Documentation

### Code Documentation
- Comprehensive docstrings for all functions
- Type hints for all parameters and returns
- Scientific citations in code comments
- Architecture and design documentation

### User Guides
- Command-line usage examples
- Error resolution guides
- Best practices for data preparation
- Troubleshooting common issues

## 🤝 Contributing

### Development Standards
- Follow established project structure
- Implement comprehensive error handling
- Include scientific citations for all constants
- Write tests for all new functionality
- Maintain data validation and cleanup

### Code Quality
- No placeholder implementations
- Complete error handling paths
- Proper resource management
- Comprehensive logging
- Performance optimization

## 📄 License

This project follows scientific software development standards with:
- Open source availability
- Scientific reproducibility requirements
- Proper attribution and citation
- Community contribution guidelines

## 🆘 Support

### Getting Help
- Check error logs for detailed information
- Review scientific citations in error messages
- Consult Ensembl API documentation
- Verify configuration parameters

### Common Issues
- **API Connection**: Ensure internet connectivity
- **File Format**: Verify CSV file format and encoding
- **Column Names**: Check column naming conventions
- **Coordinate Data**: Validate genomic coordinates

### Troubleshooting
1. **Test API Connection**: Use `--test-api` flag
2. **Check File Format**: Ensure CSV files are properly formatted
3. **Verify Columns**: Confirm required columns are present
4. **Review Logs**: Check log files for detailed error information

---

**Note**: This application is designed for scientific research and requires proper understanding of genomic analysis concepts. Always validate results against established scientific databases and literature.

## 📊 Example Usage Scenarios

### Scenario 1: Single BLAST Results File
```bash
# Process a single BLAST results file
python gene_annotation_script.py --file blast_results.csv

# Output: blast_results_annotated.csv
```

### Scenario 2: Multiple Files with Custom Output
```bash
# Process multiple files and save to results directory
python gene_annotation_script.py --file file1.csv file2.csv --output results/

# Output: results/file1_annotated.csv, results/file2_annotated.csv
```

### Scenario 3: Directory Processing
```bash
# Process all CSV files in data directory
python gene_annotation_script.py --directory data/ --pattern "*.csv"

# Output: Annotated versions of all CSV files in data/
```

### Scenario 4: Debug Mode
```bash
# Process with detailed logging
python gene_annotation_script.py --file test.csv --log-level DEBUG

# Output: Detailed logs and annotated file
```
