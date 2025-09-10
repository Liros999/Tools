# Interactive Volcano App

A comprehensive web application for interactive volcano plot visualization and analysis of gene expression data, with advanced filtering, coloring, and neighborhood exploration capabilities.

## 🚀 Features

### Core Functionality
- **Interactive Volcano Plots**: Dynamic visualization of gene expression data
- **Multi-Category Filtering**: Select and visualize genes from multiple biological categories
- **Intersection Analysis**: Identify genes that appear in multiple selected categories
- **Draggable Gene Labels**: Interactive labeling system for selected genes
- **Export Capabilities**: Download plots as PNG images with custom filenames

### Advanced Features
- **Pathway-Based Coloring**: Color genes based on metabolic pathway membership
- **Neighborhood Exploration**: Explore protein-protein interaction networks with BFS-based layering
- **Regulon Analysis**: Filter genes by transcriptional regulons
- **Protein Complex Highlighting**: Visualize protein complex memberships
- **Multi-Category Intersections**: "ALL categories" logic for finding shared genes

### Technical Features
- **Real-time Updates**: Live plot updates without page refresh
- **Responsive Design**: Mobile-friendly interface
- **Performance Optimized**: Efficient data handling and rendering
- **API Integration**: SubtiWiki, UniProt, and AlphaFold API support

## 🏗️ Architecture

### Backend (Flask)
- **Flask Framework**: Modern web application backend
- **Modular Structure**: Organized into core, API, and data management modules
- **Caching System**: Local JSON caching for external API data
- **Graph Algorithms**: BFS implementation for neighborhood exploration

### Frontend (JavaScript)
- **Plotly.js**: Interactive plotting library
- **Modular JavaScript**: Organized into specialized modules
- **State Management**: Centralized state management system
- **Event Handling**: Robust event management for user interactions

### Data Sources
- **SubtiWiki API**: Bacillus subtilis biological data
- **UniProt API**: Protein sequence and annotation data
- **AlphaFold API**: Protein structure prediction data
- **Local CSV**: User-uploaded gene expression data

## 📁 Project Structure

```
InteractiveVolcanoApp/
├── app.py                          # Main Flask application entry point
├── requirements.txt                # Python dependencies
├── README.md                      # This comprehensive documentation
├── Subtiwiki_API.txt             # SubtiWiki API documentation
├── mmc4_CueR_with_Regulators.csv # Sample gene expression data
├── Ca_D1_3_6_data_with_Regulators.csv # Additional sample data
├── src/                           # Backend source code
│   ├── core/                     # Core application logic
│   ├── api/                      # API endpoints and services
│   └── data/                     # Data management and processing
├── static/                        # Frontend assets
│   ├── js/                       # JavaScript modules
│   ├── css/                      # Stylesheets
│   └── api_cache/                # Cached API data
├── templates/                     # HTML templates
├── config/                        # Configuration files
├── docs/                         # Documentation
└── logs/                         # Application logs
```

## 🛠️ Installation

### Prerequisites
- Python 3.8+
- pip package manager
- Modern web browser

### Setup
1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd InteractiveVolcanoApp
   ```

2. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the application**
   ```bash
   python app.py
   ```

4. **Access the application**
   - Open your browser and navigate to `http://localhost:5000`
   - Upload your gene expression data (CSV format)
   - Start exploring with interactive filters and visualizations

## 📊 Data Format

### Input CSV Requirements
The application expects CSV files with the following columns:
- **Gene names**: Gene identifiers or symbols
- **Log2 Fold Change**: Expression fold change values
- **P-values**: Statistical significance values
- **Optional**: Additional metadata columns

### Sample Data
Two sample datasets are included:
- `mmc4_CueR_with_Regulators.csv`: CueR regulator analysis
- `Ca_D1_3_6_data_with_Regulators.csv`: Calcium-dependent analysis

## 🔧 Configuration

### Environment Variables
- `FLASK_ENV`: Set to 'development' for debug mode
- `FLASK_DEBUG`: Enable/disable debug features

### API Configuration
- **SubtiWiki**: Configured for Bacillus subtilis data
- **UniProt**: Protein annotation retrieval
- **AlphaFold**: Structure prediction integration

## 🎯 Usage Guide

### Basic Workflow
1. **Upload Data**: Select and upload your CSV file
2. **Configure Plot**: Set plot title and parameters
3. **Apply Filters**: Select biological categories of interest
4. **Explore Results**: Use interactive features to analyze data
5. **Export Results**: Download plots and data as needed

### Advanced Features
- **Neighborhood Exploration**: Enter a gene name and adjust layer depth
- **Pathway Coloring**: Switch to pathway-based visualization mode
- **Multi-Category Analysis**: Select multiple categories to find intersections
- **Label Management**: Toggle and drag gene labels for better visualization

## 🔍 API Endpoints

### Core Endpoints
- `POST /upload`: File upload and processing
- `GET /analysis`: Main analysis interface
- `GET /api/selection/genes`: Gene selection data
- `GET /api/neighborhood/bfs`: Neighborhood exploration
- `GET /api/pathways/list`: Metabolic pathway data
- `GET /api/regulons/list`: Transcriptional regulon data

### External APIs
- **SubtiWiki**: Biological data for Bacillus subtilis
- **UniProt**: Protein sequence and annotation data
- **AlphaFold**: Protein structure predictions

## 🧪 Testing

### Backend Testing
- Unit tests for core functionality
- API endpoint testing
- Data processing validation

### Frontend Testing
- Interactive plot functionality
- User interface responsiveness
- Cross-browser compatibility

## 🚨 Troubleshooting

### Common Issues
1. **Plot not loading**: Check browser console for JavaScript errors
2. **API failures**: Verify network connectivity and API endpoints
3. **Data not displaying**: Ensure CSV format matches requirements
4. **Performance issues**: Check browser performance and data size

### Debug Mode
Enable debug mode for detailed error messages and logging:
```python
app.run(debug=True)
```

## 📈 Performance

### Optimization Features
- **Lazy Loading**: Load data only when needed
- **Caching**: Local storage of frequently accessed data
- **Efficient Rendering**: Optimized Plotly.js configurations
- **Memory Management**: Proper cleanup of large datasets

### Scalability
- **Modular Architecture**: Easy to extend and modify
- **API Design**: RESTful endpoints for external integration
- **Data Handling**: Efficient processing of large datasets

## 🤝 Contributing

### Development Guidelines
1. **Code Style**: Follow PEP 8 Python guidelines
2. **Documentation**: Update README and inline comments
3. **Testing**: Add tests for new features
4. **Version Control**: Use descriptive commit messages

### Feature Requests
- Submit issues for bugs or feature requests
- Provide detailed descriptions of desired functionality
- Include sample data when possible

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- **SubtiWiki Team**: For providing comprehensive Bacillus subtilis data
- **Plotly.js Community**: For excellent plotting library
- **Flask Community**: For robust web framework
- **Scientific Community**: For biological data standards and APIs

## 📞 Support

For questions, issues, or contributions:
- **Issues**: Use the GitHub issue tracker
- **Documentation**: Check this README and inline code comments
- **Community**: Engage with the development community

---

**Last Updated**: December 2024  
**Version**: 1.0.0  
**Maintainer**: Interactive Volcano App Development Team
