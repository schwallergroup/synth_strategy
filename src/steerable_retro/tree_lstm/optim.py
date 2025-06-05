""" Module containing an objective class for Optuna optimization of classification model """
import torch
import optuna
from optuna.integration import PyTorchLightningPruningCallback
from pytorch_lightning import Trainer
from pytorch_lightning.callbacks import EarlyStopping, ModelCheckpoint
from route_distances.lstm.data import TreeClassificationDataModule  # Updated import
from route_distances.lstm.models import RouteClassificationModel  # Updated import

EPOCHS = 50
PATIENCE = 10


class ClassificationOptunaObjective:
    """
    Representation of an objective function for Optuna hyperparameter optimization
    
    :param data_path: the path to a pickle file with preprocessed trees and labels
    :param num_classes: number of classes for classification
    :param train_split: fraction of data for training
    :param val_split: fraction of data for validation
    :param test_split: fraction of data for testing
    """
    
    def __init__(
        self, 
        data_path: str, 
        num_classes: int,
        train_split: float = 0.7,
        val_split: float = 0.15,
        test_split: float = 0.15
    ) -> None:
        self._data_path = data_path
        self._num_classes = num_classes
        self._train_split = train_split
        self._val_split = val_split  
        self._test_split = test_split
        
    def __call__(self, trial: optuna.trial.Trial) -> float:
        # Suggest hyperparameters
        batch_size = trial.suggest_int("batch_size", 16, 128, step=16)
        lstm_size = trial.suggest_categorical("lstm_size", [256, 512, 1024, 2048])
        dropout_prob = trial.suggest_float("dropout", 0.1, 0.6, step=0.1)
        learning_rate = trial.suggest_float("lr", 1e-4, 1e-2, log=True)
        weight_decay = trial.suggest_float("weight_decay", 1e-5, 1e-2, log=True)
        threshold = trial.suggest_float("threshold", 0.3, 0.7, step=0.05)
        
        # Create data module
        data_module = TreeClassificationDataModule(
            data_path=self._data_path,
            num_classes=self._num_classes,
            batch_size=batch_size,
            train_split=self._train_split,
            val_split=self._val_split,
            test_split=self._test_split,
        )
        
        # Create model with hyperparameters
        model = RouteClassificationModel(
            num_classes=self._num_classes,
            lstm_size=lstm_size,
            dropout_prob=dropout_prob,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            threshold=threshold,
        )
        
        # Setup callbacks
        pruning_callback = PyTorchLightningPruningCallback(
            trial, monitor="val_monitor"
        )
        
        early_stopping = EarlyStopping(
            monitor="val_monitor",
            patience=PATIENCE,
            mode="max",  # We want to maximize F1 score
            verbose=False
        )
        
        checkpoint_callback = ModelCheckpoint(
            monitor="val_monitor",
            mode="max",
            save_top_k=1,
            verbose=False
        )
        
        # Setup trainer
        trainer = Trainer(
            accelerator="auto",
            devices="auto",
            logger=False,  # Disable logging for hyperparameter search
            callbacks=[pruning_callback, early_stopping, checkpoint_callback],
            max_epochs=EPOCHS,
            enable_progress_bar=False,  # Reduce output during search
            enable_model_summary=False,
        )
        
        # Train model
        try:
            trainer.fit(model, datamodule=data_module)
            
            # Return the best validation metric
            if trainer.callback_metrics:
                return trainer.callback_metrics["val_monitor"].item()
            else:
                # Fallback if no metrics available
                return 0.0
                
        except Exception as e:
            # Handle any training failures
            print(f"Trial failed with error: {e}")
            return 0.0


def optimize_hyperparameters(
    data_path: str,
    num_classes: int,
    n_trials: int = 100,
    study_name: str = "route_classification",
    direction: str = "maximize",
    train_split: float = 0.7,
    val_split: float = 0.15,
    test_split: float = 0.15
) -> optuna.Study:
    """
    Run hyperparameter optimization
    
    :param data_path: Path to the data pickle file
    :param num_classes: Number of classes
    :param n_trials: Number of optimization trials
    :param study_name: Name of the study
    :param direction: Optimization direction ("maximize" or "minimize")
    :param train_split: Training data fraction
    :param val_split: Validation data fraction  
    :param test_split: Test data fraction
    :return: Completed optuna study
    """
    
    # Create objective
    objective = ClassificationOptunaObjective(
        data_path=data_path,
        num_classes=num_classes,
        train_split=train_split,
        val_split=val_split,
        test_split=test_split
    )
    
    # Create study
    study = optuna.create_study(
        direction=direction,
        study_name=study_name,
        pruner=optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=10)
    )
    
    # Optimize
    study.optimize(objective, n_trials=n_trials)
    
    # Print results
    print("Study statistics: ")
    print(f"  Number of finished trials: {len(study.trials)}")
    print(f"  Number of pruned trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.PRUNED])}")
    print(f"  Number of complete trials: {len([t for t in study.trials if t.state == optuna.trial.TrialState.COMPLETE])}")
    
    print("Best trial:")
    trial = study.best_trial
    print(f"  Value: {trial.value}")
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")
    
    return study


# Example usage
if __name__ == "__main__":
    # Example of how to use the optimization
    study = optimize_hyperparameters(
        data_path="path/to/your/classification_data.pkl",
        num_classes=10,  # Adjust based on your problem
        n_trials=50,
    )
    
    # Save study results
    import joblib
    joblib.dump(study, f"optuna_study_{study.study_name}.pkl")