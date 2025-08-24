"""
Core genome analysis classes and functionality.
Contains the main GenomeAnalyzer class and related core components.
"""

from .genome_analyzer import GenomeAnalyzer
from .resource_manager import ResourceManager
from .exceptions import GenomeAnalysisError

__all__ = [
    'GenomeAnalyzer',
    'ResourceManager', 
    'GenomeAnalysisError'
]
