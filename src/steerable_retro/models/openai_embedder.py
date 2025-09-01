# src/steerable_retro/models/openai_embedder.py

from typing import List
import numpy as np
from openai import OpenAI
from .base_embedder import Embedder  # Direct import is clearer

class OpenAIEmbedder(Embedder):
    """Embedder for OpenAI models like text-embedding-ada-002, 3-small, 3-large."""
    def __init__(self, model_name: str, api_key: str = None):
        super().__init__(model_name)
        self.client = OpenAI(api_key=api_key)
        
        # Store dimension in a private variable
        if "3-large" in model_name:
            self._embedding_dim = 3072
        else:
            self._embedding_dim = 1536

    @property
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings."""
        return self._embedding_dim

    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """Generates embeddings for a list of texts in a batch."""
        if not texts:
            return np.array([])
            
        texts = [t.replace("\n", " ") for t in texts]
        try:
            # The 'dimensions' parameter is only for v3 models
            if self.model_name == "text-embedding-3-large":
                 response = self.client.embeddings.create(input=texts, model=self.model_name, dimensions=self._embedding_dim)
            else:
                 response = self.client.embeddings.create(input=texts, model=self.model_name)
            embeddings = [item.embedding for item in response.data]
            return np.array(embeddings)
        except Exception as e:
            print(f"Error getting OpenAI embeddings: {e}")
            return np.zeros((len(texts), self.embedding_dim))