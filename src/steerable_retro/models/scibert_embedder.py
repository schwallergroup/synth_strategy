# src/steerable_retro/models/scibert_embedder.py

from typing import List
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
from .base_embedder import Embedder # Direct import is clearer

class SciBertEmbedder(Embedder):
    """Embedder for SciBERT models using the transformers library directly."""
    def __init__(self, model_name: str = 'allenai/scibert_scivocab_uncased'):
        super().__init__(model_name)
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"SciBERT using device: {self.device}")
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = AutoModel.from_pretrained(self.model_name).to(self.device)
        self._embedding_dim = self.model.config.hidden_size

    @property
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings."""
        return self._embedding_dim

    def _mean_pooling(self, model_output, attention_mask):
        token_embeddings = model_output[0]
        input_mask_expanded = attention_mask.unsqueeze(-1).expand(token_embeddings.size()).float()
        sum_embeddings = torch.sum(token_embeddings * input_mask_expanded, 1)
        sum_mask = torch.clamp(input_mask_expanded.sum(1), min=1e-9)
        return sum_embeddings / sum_mask

    @torch.no_grad()
    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """Generates embeddings for a list of texts in a batch."""
        if not texts:
            return np.array([])
            
        encoded_input = self.tokenizer(
            texts, padding=True, truncation=True, return_tensors='pt'
        ).to(self.device)
        model_output = self.model(**encoded_input)
        sentence_embeddings = self._mean_pooling(model_output, encoded_input['attention_mask'])
        return sentence_embeddings.cpu().numpy()