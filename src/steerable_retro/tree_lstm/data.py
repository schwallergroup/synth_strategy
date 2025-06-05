""" Module containing classes for loading and generating data for classification model training """
import pickle
import random
import multiprocessing
from typing import List, Tuple, Set, Union, Dict, Any
import numpy as np
import torch
from pytorch_lightning import LightningDataModule
from torch.utils.data import Dataset, DataLoader

import route_distances.lstm.defaults as defaults
from route_distances.lstm.utils import collate_trees


class ClassificationTreeDataset(Dataset):
    """Dataset for tree-based multiclass classification"""

    def __init__(self, trees: List[Dict], labels: List[List[int]], num_classes: int):
        """
        Initialize the dataset
        
        :param trees: List of preprocessed trees
        :param labels: List of label lists (each tree can have 1-5 classes)
        :param num_classes: Total number of classes
        """
        assert len(trees) == len(labels), "Number of trees must match number of label sets"
        self.trees = trees
        self.labels = labels
        self.num_classes = num_classes

    def __len__(self):
        return len(self.trees)

    def __getitem__(self, idx):
        tree = self.trees[idx]
        label_list = self.labels[idx]
        
        # Convert label list to one-hot encoding
        one_hot = np.zeros(self.num_classes, dtype=np.float32)
        for label in label_list:
            if 0 <= label < self.num_classes:
                one_hot[label] = 1.0
                
        return {
            "tree": tree,
            "labels": one_hot
        }


def collate_classification_batch(batch: List[Dict]) -> Dict:
    """
    Collate function for classification batches
    
    :param batch: List of samples from the dataset
    :return: Collated batch
    """
    trees = [sample["tree"] for sample in batch]
    labels = torch.stack([torch.tensor(sample["labels"]) for sample in batch])
    
    collated_trees = collate_trees(trees)
    
    return {
        "tree": collated_trees,
        "labels": labels
    }


class TreeClassificationDataModule(LightningDataModule):
    """PyTorch Lightning datamodule for tree classification"""

    def __init__(
        self,
        data_path: str,  # Path to pickle file with trees and labels
        num_classes: int,
        batch_size: int = defaults.BATCH_SIZE,
        train_split: float = 0.7,
        val_split: float = 0.15,
        test_split: float = 0.15,
        random_seed: int = defaults.SPLIT_SEED,
        shuffle: bool = True,
    ) -> None:
        super().__init__()
        
        assert abs(train_split + val_split + test_split - 1.0) < 1e-6, \
            "Train, validation, and test splits must sum to 1.0"
            
        self._data_path = data_path
        self._num_classes = num_classes
        self._batch_size = batch_size
        self._train_split = train_split
        self._val_split = val_split
        self._test_split = test_split
        self._random_seed = random_seed
        self._shuffle = shuffle
        self._num_workers = min(multiprocessing.cpu_count(), 8)

        self.train_dataset = None
        self.val_dataset = None
        self.test_dataset = None

    def setup(self, stage: str = None) -> None:
        """Load and split the data"""
        with open(self._data_path, "rb") as f:
            data = pickle.load(f)
        
        # Expected format: {"trees": [...], "labels": [...]}
        trees = data["trees"]
        labels = data["labels"]
        
        assert len(trees) == len(labels), "Mismatch between number of trees and labels"
        
        # Create indices and shuffle
        indices = list(range(len(trees)))
        random.seed(self._random_seed)
        random.shuffle(indices)
        
        # Calculate split points
        n_total = len(indices)
        n_train = int(n_total * self._train_split)
        n_val = int(n_total * self._val_split)
        
        train_indices = indices[:n_train]
        val_indices = indices[n_train:n_train + n_val]
        test_indices = indices[n_train + n_val:]
        
        # Create datasets
        self.train_dataset = ClassificationTreeDataset(
            [trees[i] for i in train_indices],
            [labels[i] for i in train_indices],
            self._num_classes
        )
        
        self.val_dataset = ClassificationTreeDataset(
            [trees[i] for i in val_indices],
            [labels[i] for i in val_indices],
            self._num_classes
        )
        
        self.test_dataset = ClassificationTreeDataset(
            [trees[i] for i in test_indices],
            [labels[i] for i in test_indices],
            self._num_classes
        )
        
        print("=== Data split ===")
        print(f"Training dataset: {len(self.train_dataset)} samples")
        print(f"Validation dataset: {len(self.val_dataset)} samples")
        print(f"Test dataset: {len(self.test_dataset)} samples")
        print(f"Total classes: {self._num_classes}")

    def train_dataloader(self) -> DataLoader:
        return DataLoader(
            self.train_dataset,
            batch_size=self._batch_size,
            collate_fn=collate_classification_batch,
            shuffle=self._shuffle,
            num_workers=self._num_workers,
            pin_memory=True
        )

    def val_dataloader(self) -> DataLoader:
        return DataLoader(
            self.val_dataset,
            batch_size=self._batch_size,
            collate_fn=collate_classification_batch,
            shuffle=False,
            num_workers=self._num_workers,
            pin_memory=True
        )

    def test_dataloader(self) -> DataLoader:
        return DataLoader(
            self.test_dataset,
            batch_size=self._batch_size,
            collate_fn=collate_classification_batch,
            shuffle=False,
            num_workers=self._num_workers,
            pin_memory=True
        )

    def predict_dataloader(self) -> DataLoader:
        """Dataloader for predictions (uses test set)"""
        return self.test_dataloader()


# Utility function to create the expected data format
def create_classification_data_pickle(
    trees: List[Dict],
    labels: List[List[int]],
    output_path: str
) -> None:
    """
    Create a pickle file with the expected format for the classification data module
    
    :param trees: List of preprocessed trees (output from preprocess_reaction_tree)
    :param labels: List of label lists, where each inner list contains 1-5 class indices
    :param output_path: Path where to save the pickle file
    """
    data = {
        "trees": trees,
        "labels": labels
    }
    
    with open(output_path, "wb") as f:
        pickle.dump(data, f)
    
    print(f"Saved {len(trees)} trees with labels to {output_path}")