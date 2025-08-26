"""
B.sub_Analyzer.src package
Contains the core functionality for B. subtilis gene analysis.
"""

from .B_Sub2_0 import BSubAnalyzer
from .gene_coordinates import GeneData

__all__ = ['BSubAnalyzer', 'GeneData'] 