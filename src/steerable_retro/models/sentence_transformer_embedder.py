# src/steerable_retro/models/sentence_transformer_embedder.py

from typing import List
import numpy as np
from .base_embedder import Embedder # Direct import is clearer

# pip install sentence-transformers
from sentence_transformers import SentenceTransformer

class SentenceTransformerEmbedder(Embedder):
    """Embedder for models loaded via the sentence-transformers library."""
    def __init__(self, model_name: str):
        super().__init__(model_name)
        print(f"Loading sentence-transformer model: {model_name}...")
        self.model = SentenceTransformer(model_name)
        self._embedding_dim = self.model.get_sentence_embedding_dimension()
        print("Model loaded.")

    @property
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings."""
        return self._embedding_dim

    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """Generates embeddings for a list of texts in a batch."""
        if not texts:
            return np.array([])

        # The library handles batching and device placement automatically
        embeddings = self.model.encode(
            texts,
            convert_to_numpy=True,
            show_progress_bar=False  # Good for non-interactive scripts
        )
        return np.array(embeddings)