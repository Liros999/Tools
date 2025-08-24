"""
Custom exceptions for genome analysis operations.
Provides specific error types for different failure scenarios.
"""


class GenomeAnalysisError(Exception):
    """Custom exception for genome analysis errors."""
    pass


# Export exception
__all__ = ['GenomeAnalysisError']
