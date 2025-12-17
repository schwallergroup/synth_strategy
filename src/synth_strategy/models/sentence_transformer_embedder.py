from __future__ import annotations

from typing import List, Optional

import numpy as np
import torch
from transformers import AutoModel, AutoTokenizer

from .base_embedder import Embedder

class SentenceTransformerEmbedder(Embedder):
    """
    Inference-only embedder for SentenceTransformer-style models.

    This intentionally uses `transformers` directly (instead of importing the
    `sentence_transformers` package), because recent `sentence_transformers`
    releases pull in training-time dependencies (e.g. `accelerate`, `datasets`,
    `pyarrow`) at import time, which can break lightweight environments.
    """

    def __init__(
        self,
        model_name: str,
        *,
        device: Optional[str] = None,
        batch_size: int = 32,
        max_length: int = 256,
        local_files_only: bool = False,
    ):
        super().__init__(model_name)
        model_id = model_name if "/" in model_name else f"sentence-transformers/{model_name}"

        self.device = torch.device(device or ("cuda" if torch.cuda.is_available() else "cpu"))
        self.batch_size = int(batch_size)
        self.max_length = int(max_length)

        print(f"Loading transformer model: {model_id} (device={self.device})...")
        self.tokenizer = AutoTokenizer.from_pretrained(model_id, local_files_only=local_files_only)
        self.model = AutoModel.from_pretrained(model_id, local_files_only=local_files_only)
        self.model.eval()
        self.model.to(self.device)
        self._embedding_dim = int(getattr(self.model.config, "hidden_size", 0))
        print("Model loaded.")

    @property
    def embedding_dim(self) -> int:
        """Returns the dimension of the embeddings."""
        return self._embedding_dim

    def get_embeddings(self, texts: List[str]) -> np.ndarray:
        """Generates embeddings for a list of texts in a batch."""
        if not texts:
            return np.array([])

        all_embeddings: list[np.ndarray] = []
        with torch.no_grad():
            for start in range(0, len(texts), self.batch_size):
                batch_texts = texts[start : start + self.batch_size]
                inputs = self.tokenizer(
                    batch_texts,
                    padding=True,
                    truncation=True,
                    max_length=self.max_length,
                    return_tensors="pt",
                )
                inputs = {k: v.to(self.device) for k, v in inputs.items()}

                model_output = self.model(**inputs)
                token_embeddings = model_output.last_hidden_state
                attention_mask = inputs["attention_mask"].unsqueeze(-1).expand(token_embeddings.size()).float()

                summed = (token_embeddings * attention_mask).sum(dim=1)
                counts = attention_mask.sum(dim=1).clamp(min=1e-9)
                mean_pooled = summed / counts

                normalized = torch.nn.functional.normalize(mean_pooled, p=2, dim=1)
                all_embeddings.append(normalized.cpu().numpy())

        return np.vstack(all_embeddings) if all_embeddings else np.array([])
