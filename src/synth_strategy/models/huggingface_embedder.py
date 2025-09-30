# src/steerable_retro/models/scibert_embedder.py
from typing import List
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
from .base_embedder import Embedder # Direct import is clearer

class HuggingfaceEmbedder(Embedder):
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
        """Generates embeddings for a list of texts iteratively to minimize GPU memory usage."""
        if not texts:
            return np.array([])
        
        embeddings_list = []
        
        for text in texts:
            # Process one text at a time
            encoded_input = self.tokenizer(
                [text], padding=True, truncation=True, return_tensors='pt'
            ).to(self.device)
            
            # Get embedding for single text
            model_output = self.model(**encoded_input)
            sentence_embedding = self._mean_pooling(model_output, encoded_input['attention_mask'])
            
            # Immediately move to CPU and convert to numpy
            embedding_cpu = sentence_embedding.cpu().numpy()
            embeddings_list.append(embedding_cpu)
            
            # Clear GPU cache for this iteration
            del encoded_input, model_output, sentence_embedding
            if self.device == "cuda":
                torch.cuda.empty_cache()
        
        # Concatenate all embeddings
        return np.vstack(embeddings_list)