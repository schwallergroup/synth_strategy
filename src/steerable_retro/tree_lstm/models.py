""" Module containing the LSTM-based model for multiclass classification of routes """
from typing import List, Tuple, Dict, Any
import torch
import pytorch_lightning as lightning
from treelstm import TreeLSTM as TreeLSTMBase
from torchmetrics import Accuracy, F1Score, AUROC
import torch.nn.functional as F

import route_distances.lstm.defaults as defaults
from route_distances.lstm.utils import accumulate_stats
from route_distances.utils.type_utils import StrDict


class _TreeLstmWithPreCompression(torch.nn.Module):
    def __init__(self, fp_size: int, lstm_size: int, dropout_prob: float) -> None:
        super().__init__()
        self._compression = torch.nn.Sequential(
            torch.nn.Linear(fp_size, lstm_size),
            torch.nn.ReLU(),
            torch.nn.Dropout(p=dropout_prob),
            torch.nn.Linear(lstm_size, lstm_size),
            torch.nn.ReLU(),
        )
        self._tree_lstm = TreeLSTMBase(lstm_size, lstm_size)

    def forward(self, tree_batch: StrDict) -> torch.Tensor:
        """
        Forward pass

        :param tree_batch: collated trees from the `route_distances.utils.collate_trees` function.
        :return: the LSTM representation of the first nodes
        """
        features = self._compression(tree_batch["features"])
        lstm_output, _ = self._tree_lstm(
            features,
            tree_batch["node_order"],
            tree_batch["adjacency_list"],
            tree_batch["edge_order"],
        )
        # Only save value of top-node
        lstm_output = torch.stack(
            [t[0, :] for t in torch.split(lstm_output, tree_batch["tree_sizes"], dim=0)]
        )
        return lstm_output


class RouteClassificationModel(lightning.LightningModule):
    """
    Model for multiclass classification of synthesis routes

    :param fp_size: the length of the fingerprint vector
    :param lstm_size: the size of the LSTM cell
    :param num_classes: the number of classes for classification
    :param dropout_prob: the dropout probability
    :param learning_rate: the initial learning rate of the optimizer
    :param weight_decay: weight decay factor of the optimizer
    :param threshold: threshold for binary classification in multilabel setting
    """

    def __init__(
        self,
        num_classes: int,
        fp_size: int = defaults.FP_SIZE,
        lstm_size: int = defaults.LSTM_SIZE,
        dropout_prob: float = defaults.DROPOUT_PROB,
        learning_rate: float = defaults.LEARNING_RATE,
        weight_decay: float = defaults.WEIGHT_DECAY,
        threshold: float = 0.5,
    ) -> None:
        super().__init__()
        self.save_hyperparameters()
        
        self._tree_lstm = _TreeLstmWithPreCompression(fp_size, lstm_size, dropout_prob)
        
        # Classification head
        self._classifier = torch.nn.Sequential(
            torch.nn.Linear(lstm_size, lstm_size // 2),
            torch.nn.ReLU(),
            torch.nn.Dropout(p=dropout_prob),
            torch.nn.Linear(lstm_size // 2, num_classes)
        )
        
        self._loss_func = torch.nn.BCEWithLogitsLoss()  # For multilabel classification
        self._threshold = threshold
        self._num_classes = num_classes
        
        # Metrics
        self._accuracy = Accuracy(task="multilabel", num_labels=num_classes, threshold=threshold)
        self._f1_micro = F1Score(task="multilabel", num_labels=num_classes, average="micro")
        self._f1_macro = F1Score(task="multilabel", num_labels=num_classes, average="macro")
        self._auroc = AUROC(task="multilabel", num_labels=num_classes)
        
        self._lr = learning_rate
        self._weight_decay = weight_decay

    def forward(self, tree_data: StrDict) -> torch.Tensor:
        """
        Forward pass through the model

        :param tree_data: collated trees from the `route_distances.utils.collate_trees` function.
        :return: the class logits
        """
        lstm_enc = self._tree_lstm(tree_data)
        logits = self._classifier(lstm_enc)
        return logits

    def training_step(self, batch: StrDict, batch_idx: int) -> torch.Tensor:
        """
        One step in the training loop

        :param batch: batch containing tree data and labels
        :param batch_idx: batch index
        :return: the loss tensor
        """
        logits = self(batch["tree"])
        loss = self._loss_func(logits, batch["labels"].float())
        
        # Calculate probabilities for metrics
        probs = torch.sigmoid(logits)
        
        # Log metrics
        self.log("train_loss", loss, prog_bar=True)
        self.log("train_acc", self._accuracy(probs, batch["labels"]), prog_bar=True)
        self.log("train_f1_micro", self._f1_micro(probs, batch["labels"]))
        self.log("train_f1_macro", self._f1_macro(probs, batch["labels"]))
        
        return loss

    def validation_step(self, batch: StrDict, batch_idx: int) -> Dict[str, torch.Tensor]:
        """
        One step in the validation loop

        :param batch: batch containing tree data and labels
        :param batch_idx: batch index
        :return: the validation metrics
        """
        return self._val_and_test_step(batch, "val")

    def validation_epoch_end(self, outputs: List[Dict[str, torch.Tensor]]) -> None:
        """Log the average validation metrics"""
        self._log_average_metrics(outputs, "val")

    def test_step(self, batch: StrDict, batch_idx: int) -> Dict[str, torch.Tensor]:
        """
        One step in the test loop

        :param batch: batch containing tree data and labels
        :param batch_idx: batch index
        :return: the test metrics
        """
        return self._val_and_test_step(batch, "test")

    def test_epoch_end(self, outputs: List[Dict[str, torch.Tensor]]) -> None:
        """Log the average test metrics"""
        self._log_average_metrics(outputs, "test")

    def configure_optimizers(self) -> Tuple[List[torch.optim.Adam], List[StrDict]]:
        """Setup the Adam optimiser and scheduler"""
        optim = torch.optim.Adam(
            self.parameters(), lr=self._lr, weight_decay=self._weight_decay
        )
        scheduler = {
            "scheduler": torch.optim.lr_scheduler.ReduceLROnPlateau(optim, patience=5),
            "monitor": "val_loss",
        }
        return [optim], [scheduler]

    def predict_step(self, batch: StrDict, batch_idx: int) -> torch.Tensor:
        """
        Prediction step

        :param batch: batch containing tree data
        :param batch_idx: batch index
        :return: predicted probabilities
        """
        logits = self(batch["tree"])
        probs = torch.sigmoid(logits)
        return probs

    def _val_and_test_step(self, batch: StrDict, prefix: str) -> Dict[str, torch.Tensor]:
        """Common validation and test step logic"""
        logits = self(batch["tree"])
        loss = self._loss_func(logits, batch["labels"].float())
        probs = torch.sigmoid(logits)
        
        return {
            f"{prefix}_loss": loss,
            f"{prefix}_acc": self._accuracy(probs, batch["labels"]),
            f"{prefix}_f1_micro": self._f1_micro(probs, batch["labels"]),
            f"{prefix}_f1_macro": self._f1_macro(probs, batch["labels"]),
            f"{prefix}_auroc": self._auroc(probs, batch["labels"]),
        }

    def _log_average_metrics(self, outputs: List[Dict[str, torch.Tensor]], prefix: str) -> None:
        """Log average metrics across epoch"""
        avg_loss = torch.stack([x[f"{prefix}_loss"] for x in outputs]).mean()
        avg_acc = torch.stack([x[f"{prefix}_acc"] for x in outputs]).mean()
        avg_f1_micro = torch.stack([x[f"{prefix}_f1_micro"] for x in outputs]).mean()
        avg_f1_macro = torch.stack([x[f"{prefix}_f1_macro"] for x in outputs]).mean()
        avg_auroc = torch.stack([x[f"{prefix}_auroc"] for x in outputs]).mean()
        
        self.log(f"{prefix}_loss", avg_loss, prog_bar=True)
        self.log(f"{prefix}_acc", avg_acc, prog_bar=True)
        self.log(f"{prefix}_f1_micro", avg_f1_micro)
        self.log(f"{prefix}_f1_macro", avg_f1_macro)
        self.log(f"{prefix}_auroc", avg_auroc)
        
        # For monitoring/early stopping
        self.log(f"{prefix}_monitor", avg_f1_macro)