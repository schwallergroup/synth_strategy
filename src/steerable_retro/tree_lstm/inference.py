""" Module containing class to make classification predictions """
import numpy as np
import torch
from typing import List, Dict, Tuple
from route_distances.lstm.features import preprocess_reaction_tree
from route_distances.lstm.utils import collate_trees
from route_distances.lstm.models import RouteClassificationModel  # Updated import
from route_distances.utils.type_utils import RouteList


class ClassificationInferenceHelper:
    """
    Helper class for making classification predictions using LSTM model
    
    :param model_path: the path to the model checkpoint file
    :param threshold: threshold for binary classification (default: 0.5)
    """
    
    def __init__(self, model_path: str, threshold: float = 0.5) -> None:
        self._model = RouteClassificationModel.load_from_checkpoint(model_path)
        self._model.eval()
        self._threshold = threshold
        
    def predict_probabilities(self, routes: RouteList) -> np.ndarray:
        """
        Predict class probabilities for a list of routes
        
        :param routes: List of routes in dictionary format
        :return: numpy array of shape (n_routes, n_classes) with probabilities
        """
        trees = [
            preprocess_reaction_tree(route, self._model.hparams.fp_size)
            for route in routes
        ]
        tree_data = collate_trees(trees)
        
        with torch.no_grad():
            logits = self._model(tree_data)
            probabilities = torch.sigmoid(logits)
            
        return probabilities.detach().numpy()
    
    def predict_classes(self, routes: RouteList, return_probabilities: bool = False) -> Dict:
        """
        Predict classes for a list of routes
        
        :param routes: List of routes in dictionary format
        :param return_probabilities: Whether to also return probabilities
        :return: Dictionary with predictions and optionally probabilities
        """
        probabilities = self.predict_probabilities(routes)
        
        # Apply threshold to get binary predictions
        binary_predictions = (probabilities >= self._threshold).astype(int)
        
        # Convert to list of class indices for each route
        class_predictions = []
        for i, pred in enumerate(binary_predictions):
            classes = np.where(pred == 1)[0].tolist()
            class_predictions.append(classes)
        
        result = {
            "predicted_classes": class_predictions,
            "binary_predictions": binary_predictions
        }
        
        if return_probabilities:
            result["probabilities"] = probabilities
            
        return result
    
    def predict_top_k_classes(self, routes: RouteList, k: int = 3) -> List[List[Tuple[int, float]]]:
        """
        Predict top-k classes for each route
        
        :param routes: List of routes in dictionary format
        :param k: Number of top classes to return
        :return: List of lists, where each inner list contains (class_index, probability) tuples
        """
        probabilities = self.predict_probabilities(routes)
        
        top_k_predictions = []
        for prob_row in probabilities:
            # Get indices of top-k probabilities
            top_k_indices = np.argsort(prob_row)[-k:][::-1]  # Descending order
            top_k_with_probs = [(int(idx), float(prob_row[idx])) for idx in top_k_indices]
            top_k_predictions.append(top_k_with_probs)
            
        return top_k_predictions


# Global instance cache for efficient reuse
_inst_models = {}


def get_classification_predictor(model_path: str, threshold: float = 0.5) -> ClassificationInferenceHelper:
    """
    Get a classification predictor instance (cached for efficiency)
    
    :param model_path: Path to the model checkpoint
    :param threshold: Threshold for binary classification
    :return: ClassificationInferenceHelper instance
    """
    global _inst_models
    key = f"{model_path}_{threshold}"
    
    if key not in _inst_models:
        _inst_models[key] = ClassificationInferenceHelper(model_path, threshold)
    
    return _inst_models[key]


# Convenience functions
def predict_route_classes(
    routes: RouteList, 
    model_path: str, 
    threshold: float = 0.5,
    return_probabilities: bool = False
) -> Dict:
    """
    Convenience function to predict classes for routes
    
    :param routes: List of routes in dictionary format
    :param model_path: Path to the trained model
    :param threshold: Threshold for binary classification
    :param return_probabilities: Whether to return probabilities
    :return: Dictionary with predictions
    """
    predictor = get_classification_predictor(model_path, threshold)
    return predictor.predict_classes(routes, return_probabilities)


def predict_route_probabilities(routes: RouteList, model_path: str) -> np.ndarray:
    """
    Convenience function to predict probabilities for routes
    
    :param routes: List of routes in dictionary format  
    :param model_path: Path to the trained model
    :return: Array of probabilities
    """
    predictor = get_classification_predictor(model_path)
    return predictor.predict_probabilities(routes)