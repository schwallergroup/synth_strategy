from abc import ABC, abstractmethod
from typing import List
import numpy as np

class Embedder(ABC):
    """Abstract base class for all embedder models."""
    def __init__(self, model_name: str, **kwargs):
        self.model_name = model_name

    @property
    @abstractmethod
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings."""
        pass

    @abstractmethod
    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """Generates embeddings for a list of texts."""
        pass
    
    def get_embedding(self, text: str) -> np.ndarray:
        """
        Generates an embedding for a single string of text.
        This default implementation uses the batch method.
        """
        if not text:
            return np.zeros(self.embedding_dim)
        # Call the batch method with a single-item list and get the first result
        return self.get_embeddings([text])[0]