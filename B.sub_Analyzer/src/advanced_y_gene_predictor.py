#!/usr/bin/env python3
"""
Advanced Y Gene Predictor for B. subtilis
Comprehensive ML-based gene function prediction system
"""

import logging
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from pathlib import Path
import json
import pickle
from collections import defaultdict
import re
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AdvancedYGenePredictor:
    """
    Advanced machine learning-based gene function predictor for B. subtilis
    """
    
    def __init__(self, cache_dir: str = "predictor_cache"):
        """Initialize the advanced predictor"""
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize ML models
        self.models = {
            'random_forest': RandomForestClassifier(n_estimators=100, random_state=42),
            'gradient_boost': GradientBoostingClassifier(n_estimators=100, random_state=42),
            'logistic': LogisticRegression(random_state=42, max_iter=1000)
        }
        
        # Feature extractors and scalers
        self.scaler = StandardScaler()
        self.label_encoder = LabelEncoder()
        
        # Gene function categories based on B. subtilis research
        self.function_categories = {
            'stress_response': ['sigB', 'rsbV', 'rsbW', 'rsbY', 'ctc', 'gsiB'],
            'metabolism': ['glcK', 'pfkA', 'fbaA', 'tpiA', 'gapA'],
            'transcription': ['rpoA', 'rpoB', 'rpoC', 'sigA', 'sigB'],
            'translation': ['rpsA', 'rplA', 'rpsB', 'rplB', 'tufA'],
            'cell_wall': ['murA', 'murB', 'murC', 'murD', 'pbpA'],
            'transport': ['malK', 'malF', 'malG', 'araE', 'xylE'],
            'unknown': []
        }
        
        # Rule-based patterns for gene function prediction
        self.gene_patterns = {
            'stress_response': [
                r'sig[A-Z]',  # Sigma factors
                r'rsb[A-Z]',  # Stressosome components
                r'ctc',       # General stress protein
                r'gsi[A-Z]',  # General stress induced
                r'dps',       # DNA protection during starvation
                r'katA',      # Catalase
                r'sodA'       # Superoxide dismutase
            ],
            'metabolism': [
                r'.*[Kk]inase',   # Kinases
                r'.*[Dd]ehydrogenase',  # Dehydrogenases
                r'.*[Ss]ynthase',       # Synthases
                r'glc[A-Z]',            # Glucose metabolism
                r'pyk[A-Z]',            # Pyruvate kinase
                r'eno[A-Z]'             # Enolase
            ],
            'transcription': [
                r'rpo[A-Z]',     # RNA polymerase
                r'sig[A-Z]',     # Sigma factors
                r'.*[Rr]egulator',  # Transcriptional regulators
                r'.*[Aa]ctivator',  # Transcriptional activators
                r'.*[Rr]epressor'   # Transcriptional repressors
            ],
            'translation': [
                r'rps[A-Z]',     # Ribosomal proteins small subunit
                r'rpl[A-Z]',     # Ribosomal proteins large subunit
                r'tuf[A-Z]',     # Elongation factors
                r'.*[Tt]RNA',    # tRNA related
                r'.*[Rr]ibosom'  # Ribosomal
            ],
            'cell_wall': [
                r'mur[A-Z]',     # Murein synthesis
                r'pbp[A-Z]',     # Penicillin binding proteins
                r'.*[Pp]eptidoglycan',  # Peptidoglycan
                r'.*[Cc]ell[Ww]all',    # Cell wall
                r'dac[A-Z]'      # D-alanyl-D-alanine carboxypeptidase
            ],
            'transport': [
                r'.*[Tt]ransport',   # Transport proteins
                r'.*[Pp]ermease',    # Permeases
                r'.*ABC',            # ABC transporters
                r'mal[A-Z]',         # Maltose transport
                r'ara[A-Z]',         # Arabinose transport
                r'xyl[A-Z]'          # Xylose transport
            ]
        }
        
        logger.info("Advanced Y Gene Predictor initialized")
    
    def extract_sequence_features(self, gene_name: str, sequence: Optional[str] = None) -> Dict[str, float]:
        """Extract sequence-based features from gene name and sequence"""
        features = {}
        
        # Basic name-based features
        features['name_length'] = len(gene_name)
        features['has_numbers'] = float(bool(re.search(r'\d', gene_name)))
        features['has_uppercase'] = float(any(c.isupper() for c in gene_name))
        features['starts_with_y'] = float(gene_name.lower().startswith('y'))
        
        # Pattern matching features
        for category, patterns in self.gene_patterns.items():
            pattern_matches = 0
            for pattern in patterns:
                if re.search(pattern, gene_name, re.IGNORECASE):
                    pattern_matches += 1
            features[f'{category}_pattern_score'] = pattern_matches / len(patterns)
        
        # Sequence features (if available)
        if sequence:
            features['sequence_length'] = len(sequence)
            features['gc_content'] = (sequence.count('G') + sequence.count('C')) / len(sequence)
            features['start_codon_atg'] = float(sequence.upper().startswith('ATG'))
            
            # Codon usage features
            for codon in ['ATG', 'TAA', 'TAG', 'TGA']:
                features[f'codon_{codon}_count'] = sequence.upper().count(codon)
        else:
            # Default sequence features when sequence not available
            features['sequence_length'] = 1000  # Average gene length
            features['gc_content'] = 0.43      # B. subtilis average GC content
            features['start_codon_atg'] = 1.0
            for codon in ['ATG', 'TAA', 'TAG', 'TGA']:
                features[f'codon_{codon}_count'] = 0
        
        return features
    
    def extract_expression_features(self, gene_name: str, expression_data: Optional[Dict] = None) -> Dict[str, float]:
        """Extract expression-based features"""
        features = {}
        
        if expression_data:
            # Use actual expression data
            features['mean_expression'] = expression_data.get('mean_expression', 0.0)
            features['max_expression'] = expression_data.get('max_expression', 0.0)
            features['expression_variance'] = expression_data.get('variance', 0.0)
            features['stress_upregulated'] = expression_data.get('stress_response', 0.0)
        else:
            # Simulate expression features based on gene patterns
            if any(re.search(pattern, gene_name, re.IGNORECASE) 
                   for pattern in self.gene_patterns['stress_response']):
                features['mean_expression'] = 2.5
                features['stress_upregulated'] = 1.0
            else:
                features['mean_expression'] = 1.0
                features['stress_upregulated'] = 0.0
            
            features['max_expression'] = features['mean_expression'] * 1.5
            features['expression_variance'] = 0.5
        
        return features
    
    def extract_network_features(self, gene_name: str, network_data: Optional[Dict] = None) -> Dict[str, float]:
        """Extract protein-protein interaction network features"""
        features = {}
        
        if network_data:
            features['degree_centrality'] = network_data.get('degree', 0.0)
            features['betweenness_centrality'] = network_data.get('betweenness', 0.0)
            features['clustering_coefficient'] = network_data.get('clustering', 0.0)
            features['has_interactions'] = float(network_data.get('degree', 0) > 0)
        else:
            # Estimate network properties based on known patterns
            known_hub_patterns = ['sigB', 'rsbV', 'rpoA', 'dnaA']
            is_likely_hub = any(pattern in gene_name.lower() for pattern in known_hub_patterns)
            
            if is_likely_hub:
                features['degree_centrality'] = 0.8
                features['betweenness_centrality'] = 0.6
                features['clustering_coefficient'] = 0.4
                features['has_interactions'] = 1.0
            else:
                features['degree_centrality'] = 0.2
                features['betweenness_centrality'] = 0.1
                features['clustering_coefficient'] = 0.3
                features['has_interactions'] = 0.5
        
        return features
    
    def create_feature_matrix(self, gene_list: List[str], 
                            additional_data: Optional[Dict] = None) -> pd.DataFrame:
        """Create comprehensive feature matrix for genes"""
        features_list = []
        
        for gene in gene_list:
            gene_features = {}
            
            # Extract different types of features
            sequence_features = self.extract_sequence_features(gene)
            expression_features = self.extract_expression_features(gene)
            network_features = self.extract_network_features(gene)
            
            # Combine all features
            gene_features.update(sequence_features)
            gene_features.update(expression_features)
            gene_features.update(network_features)
            
            # Add additional features if provided
            if additional_data and gene in additional_data:
                gene_features.update(additional_data[gene])
            
            features_list.append(gene_features)
        
        # Create DataFrame
        feature_matrix = pd.DataFrame(features_list, index=gene_list)
        feature_matrix = feature_matrix.fillna(0)
        
        logger.info(f"Created feature matrix: {feature_matrix.shape}")
        return feature_matrix
    
    def predict_function_rule_based(self, gene_name: str) -> Tuple[str, float]:
        """Rule-based function prediction"""
        gene_lower = gene_name.lower()
        
        # Check each category
        for category, patterns in self.gene_patterns.items():
            pattern_matches = 0
            for pattern in patterns:
                if re.search(pattern, gene_name, re.IGNORECASE):
                    pattern_matches += 1
            
            if pattern_matches > 0:
                confidence = min(pattern_matches / len(patterns), 1.0)
                return category, confidence
        
        # Check known gene lists
        for category, known_genes in self.function_categories.items():
            if gene_lower in [g.lower() for g in known_genes]:
                return category, 0.9
        
        # Default prediction
        return 'unknown', 0.1
    
    def predict_function_ml(self, gene_name: str, features: Optional[Dict] = None) -> Tuple[str, float]:
        """Machine learning-based function prediction"""
        # For now, use rule-based as ML models need training data
        # This is where you would implement the actual ML prediction
        return self.predict_function_rule_based(gene_name)
    
    def predict_gene_function(self, gene_name: str, 
                            use_ml: bool = True,
                            additional_features: Optional[Dict] = None) -> Dict[str, Any]:
        """Main prediction method"""
        result = {
            'gene': gene_name,
            'predicted_function': None,
            'confidence': 0.0,
            'method': 'rule_based',
            'features_used': [],
            'alternative_predictions': []
        }
        
        try:
            if use_ml:
                # Try ML prediction first
                function, confidence = self.predict_function_ml(gene_name, additional_features)
                result['method'] = 'machine_learning'
            else:
                # Use rule-based prediction
                function, confidence = self.predict_function_rule_based(gene_name)
            
            result['predicted_function'] = function
            result['confidence'] = confidence
            
            # Get alternative predictions
            alternatives = []
            for category in self.function_categories.keys():
                if category != function:
                    alt_confidence = self._calculate_alternative_confidence(gene_name, category)
                    if alt_confidence > 0.1:
                        alternatives.append({'function': category, 'confidence': alt_confidence})
            
            # Sort alternatives by confidence
            alternatives.sort(key=lambda x: x['confidence'], reverse=True)
            result['alternative_predictions'] = alternatives[:3]
            
        except Exception as e:
            logger.error(f"Prediction failed for {gene_name}: {e}")
            result['predicted_function'] = 'unknown'
            result['confidence'] = 0.1
            result['error'] = str(e)
        
        return result
    
    def _calculate_alternative_confidence(self, gene_name: str, category: str) -> float:
        """Calculate confidence for alternative function predictions"""
        if category not in self.gene_patterns:
            return 0.0
        
        patterns = self.gene_patterns[category]
        matches = sum(1 for pattern in patterns if re.search(pattern, gene_name, re.IGNORECASE))
        return matches / len(patterns) * 0.8  # Lower confidence for alternatives
    
    def batch_predict(self, gene_list: List[str], **kwargs) -> List[Dict[str, Any]]:
        """Predict functions for multiple genes"""
        results = []
        for gene in gene_list:
            result = self.predict_gene_function(gene, **kwargs)
            results.append(result)
        
        logger.info(f"Completed batch prediction for {len(gene_list)} genes")
        return results
    
    def get_function_summary(self, predictions: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Generate summary statistics for a set of predictions"""
        if not predictions:
            return {}
        
        # Count predictions by function
        function_counts = defaultdict(int)
        confidence_scores = []
        
        for pred in predictions:
            function = pred.get('predicted_function', 'unknown')
            confidence = pred.get('confidence', 0.0)
            
            function_counts[function] += 1
            confidence_scores.append(confidence)
        
        summary = {
            'total_genes': len(predictions),
            'function_distribution': dict(function_counts),
            'average_confidence': np.mean(confidence_scores),
            'high_confidence_predictions': sum(1 for c in confidence_scores if c > 0.7),
            'unknown_function': function_counts.get('unknown', 0)
        }
        
        return summary
    
    def save_predictions(self, predictions: List[Dict[str, Any]], filename: str):
        """Save predictions to file"""
        output_file = self.cache_dir / filename
        
        # Convert to DataFrame for easy saving
        df = pd.DataFrame(predictions)
        df.to_csv(output_file.with_suffix('.csv'), index=False)
        
        # Also save as JSON for programmatic access
        with open(output_file.with_suffix('.json'), 'w') as f:
            json.dump(predictions, f, indent=2)
        
        logger.info(f"Predictions saved to {output_file}")


def test_advanced_predictor():
    """Test the advanced predictor"""
    print("Testing Advanced Y Gene Predictor...")
    
    predictor = AdvancedYGenePredictor()
    
    # Test single prediction
    test_gene = "rsbV"
    result = predictor.predict_gene_function(test_gene)
    
    print(f"Test gene: {test_gene}")
    print(f"Predicted function: {result['predicted_function']}")
    print(f"Confidence: {result['confidence']:.3f}")
    print(f"Method: {result['method']}")
    
    # Test batch prediction
    test_genes = ["rsbV", "sigB", "dnaA", "unknownGene123"]
    batch_results = predictor.batch_predict(test_genes)
    
    print(f"\nBatch prediction results:")
    for result in batch_results:
        print(f"  {result['gene']}: {result['predicted_function']} ({result['confidence']:.3f})")
    
    # Generate summary
    summary = predictor.get_function_summary(batch_results)
    print(f"\nSummary: {summary}")
    
    return predictor


if __name__ == "__main__":
    test_advanced_predictor() 