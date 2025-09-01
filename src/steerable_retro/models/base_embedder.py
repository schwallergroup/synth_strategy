# src/steerable_retro/models/base_embedder.py

from abc import ABC, abstractmethod
from typing import List
import numpy as np

class Embedder(ABC):
    """Abstract base class for embedding models."""
    def __init__(self, model_name: str, **kwargs):
        self.model_name = model_name

    @property
    @abstractmethod
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings. Must be implemented by subclasses."""
        pass

    @abstractmethod
    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """
        Generates embeddings for a list of texts in a batch.
        This is the primary method that subclasses must implement.
        """
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