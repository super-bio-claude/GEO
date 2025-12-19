"""
Machine Learning Models for RNA-seq Classification
===================================================

This module implements various ML classifiers:
1. Logistic Regression
2. Random Forest
3. Support Vector Machine (SVM)
4. XGBoost
5. Neural Network (MLP)

With hyperparameter tuning and cross-validation.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import (
    train_test_split,
    cross_val_score,
    StratifiedKFold,
    GridSearchCV,
    RandomizedSearchCV
)
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    confusion_matrix,
    classification_report
)
from typing import Dict, List, Tuple, Optional, Any
import logging
import joblib
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class RNAseqClassifier:
    """ML classifier for RNA-seq data."""

    def __init__(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        test_size: float = 0.2,
        random_state: int = 42
    ):
        """
        Initialize classifier.

        Parameters
        ----------
        X : pd.DataFrame
            Feature matrix (samples x features)
        y : pd.Series
            Labels
        test_size : float
            Fraction of data for testing
        random_state : int
            Random seed
        """
        self.X = X
        self.y = y
        self.random_state = random_state

        # Encode labels
        self.le = LabelEncoder()
        self.y_encoded = self.le.fit_transform(y)

        # Split data
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
            X, self.y_encoded,
            test_size=test_size,
            stratify=self.y_encoded,
            random_state=random_state
        )

        logger.info(f"Train set: {len(self.X_train)} samples")
        logger.info(f"Test set: {len(self.X_test)} samples")

        self.models = {}
        self.results = {}

    def get_model_configs(self) -> Dict[str, Dict]:
        """Get default model configurations."""
        return {
            'LogisticRegression': {
                'model': LogisticRegression(
                    max_iter=1000,
                    random_state=self.random_state
                ),
                'params': {
                    'classifier__C': [0.01, 0.1, 1, 10],
                    'classifier__penalty': ['l1', 'l2'],
                    'classifier__solver': ['saga']
                }
            },
            'RandomForest': {
                'model': RandomForestClassifier(
                    random_state=self.random_state,
                    n_jobs=-1
                ),
                'params': {
                    'classifier__n_estimators': [100, 200, 500],
                    'classifier__max_depth': [5, 10, 20, None],
                    'classifier__min_samples_split': [2, 5, 10],
                    'classifier__min_samples_leaf': [1, 2, 4]
                }
            },
            'SVM': {
                'model': SVC(
                    probability=True,
                    random_state=self.random_state
                ),
                'params': {
                    'classifier__C': [0.1, 1, 10],
                    'classifier__kernel': ['linear', 'rbf'],
                    'classifier__gamma': ['scale', 'auto']
                }
            },
            'GradientBoosting': {
                'model': GradientBoostingClassifier(
                    random_state=self.random_state
                ),
                'params': {
                    'classifier__n_estimators': [100, 200],
                    'classifier__max_depth': [3, 5, 7],
                    'classifier__learning_rate': [0.01, 0.1, 0.3]
                }
            },
            'MLP': {
                'model': MLPClassifier(
                    max_iter=1000,
                    random_state=self.random_state
                ),
                'params': {
                    'classifier__hidden_layer_sizes': [(50,), (100,), (50, 50), (100, 50)],
                    'classifier__alpha': [0.0001, 0.001, 0.01],
                    'classifier__learning_rate': ['constant', 'adaptive']
                }
            }
        }

    def create_pipeline(self, model) -> Pipeline:
        """Create sklearn pipeline with scaling."""
        return Pipeline([
            ('scaler', StandardScaler()),
            ('classifier', model)
        ])

    def train_model(
        self,
        model_name: str,
        tune_hyperparameters: bool = True,
        cv_folds: int = 5,
        n_iter: int = 20
    ) -> Dict[str, Any]:
        """
        Train a single model.

        Parameters
        ----------
        model_name : str
            Name of the model to train
        tune_hyperparameters : bool
            Whether to perform hyperparameter tuning
        cv_folds : int
            Number of cross-validation folds
        n_iter : int
            Number of iterations for random search

        Returns
        -------
        Dict[str, Any]
            Training results
        """
        logger.info(f"Training {model_name}")

        configs = self.get_model_configs()
        if model_name not in configs:
            raise ValueError(f"Unknown model: {model_name}")

        config = configs[model_name]
        pipeline = self.create_pipeline(config['model'])

        if tune_hyperparameters and config['params']:
            # Use RandomizedSearchCV for efficiency
            search = RandomizedSearchCV(
                pipeline,
                config['params'],
                n_iter=min(n_iter, self._count_param_combinations(config['params'])),
                cv=StratifiedKFold(n_splits=cv_folds, shuffle=True, random_state=self.random_state),
                scoring='accuracy',
                n_jobs=-1,
                random_state=self.random_state
            )
            search.fit(self.X_train, self.y_train)
            best_model = search.best_estimator_
            best_params = search.best_params_
            cv_score = search.best_score_
        else:
            best_model = pipeline
            best_model.fit(self.X_train, self.y_train)
            best_params = {}
            cv_score = cross_val_score(
                best_model, self.X_train, self.y_train,
                cv=cv_folds, scoring='accuracy'
            ).mean()

        # Evaluate on test set
        y_pred = best_model.predict(self.X_test)
        y_prob = best_model.predict_proba(self.X_test)[:, 1] if hasattr(best_model, 'predict_proba') else None

        # Calculate metrics
        results = {
            'model_name': model_name,
            'best_params': best_params,
            'cv_score': cv_score,
            'test_accuracy': accuracy_score(self.y_test, y_pred),
            'test_precision': precision_score(self.y_test, y_pred, average='weighted'),
            'test_recall': recall_score(self.y_test, y_pred, average='weighted'),
            'test_f1': f1_score(self.y_test, y_pred, average='weighted'),
            'test_auc': roc_auc_score(self.y_test, y_prob) if y_prob is not None else None,
            'confusion_matrix': confusion_matrix(self.y_test, y_pred).tolist()
        }

        self.models[model_name] = best_model
        self.results[model_name] = results

        logger.info(f"{model_name} - Test Accuracy: {results['test_accuracy']:.4f}, "
                   f"CV Score: {cv_score:.4f}")

        return results

    def _count_param_combinations(self, params: Dict) -> int:
        """Count total parameter combinations."""
        count = 1
        for values in params.values():
            count *= len(values)
        return count

    def train_all_models(
        self,
        tune_hyperparameters: bool = True,
        cv_folds: int = 5
    ) -> pd.DataFrame:
        """
        Train all available models.

        Returns
        -------
        pd.DataFrame
            Results summary for all models
        """
        configs = self.get_model_configs()

        for model_name in configs.keys():
            try:
                self.train_model(
                    model_name,
                    tune_hyperparameters=tune_hyperparameters,
                    cv_folds=cv_folds
                )
            except Exception as e:
                logger.error(f"Error training {model_name}: {e}")

        return self.get_results_summary()

    def get_results_summary(self) -> pd.DataFrame:
        """Get summary of all model results."""
        summary = []

        for model_name, results in self.results.items():
            summary.append({
                'model': model_name,
                'cv_score': results['cv_score'],
                'test_accuracy': results['test_accuracy'],
                'test_precision': results['test_precision'],
                'test_recall': results['test_recall'],
                'test_f1': results['test_f1'],
                'test_auc': results['test_auc']
            })

        df = pd.DataFrame(summary)
        return df.sort_values('test_f1', ascending=False)

    def get_feature_importance(
        self,
        model_name: str = 'RandomForest'
    ) -> pd.DataFrame:
        """
        Get feature importance from a trained model.

        Parameters
        ----------
        model_name : str
            Model to get importance from

        Returns
        -------
        pd.DataFrame
            Feature importance scores
        """
        if model_name not in self.models:
            raise ValueError(f"Model {model_name} not trained")

        model = self.models[model_name]
        classifier = model.named_steps['classifier']

        if hasattr(classifier, 'feature_importances_'):
            importance = classifier.feature_importances_
        elif hasattr(classifier, 'coef_'):
            importance = np.abs(classifier.coef_).flatten()
        else:
            raise ValueError(f"Model {model_name} does not support feature importance")

        importance_df = pd.DataFrame({
            'feature': self.X.columns,
            'importance': importance
        }).sort_values('importance', ascending=False)

        return importance_df

    def save_model(self, model_name: str, filepath: str):
        """Save trained model to file."""
        if model_name not in self.models:
            raise ValueError(f"Model {model_name} not trained")

        joblib.dump(self.models[model_name], filepath)
        logger.info(f"Saved {model_name} to {filepath}")

    def load_model(self, model_name: str, filepath: str):
        """Load model from file."""
        self.models[model_name] = joblib.load(filepath)
        logger.info(f"Loaded {model_name} from {filepath}")

    def predict(
        self,
        X_new: pd.DataFrame,
        model_name: str = 'RandomForest'
    ) -> np.ndarray:
        """
        Make predictions on new data.

        Parameters
        ----------
        X_new : pd.DataFrame
            New samples to predict
        model_name : str
            Model to use

        Returns
        -------
        np.ndarray
            Predicted labels
        """
        if model_name not in self.models:
            raise ValueError(f"Model {model_name} not trained")

        predictions = self.models[model_name].predict(X_new)
        return self.le.inverse_transform(predictions)


def nested_cross_validation(
    X: pd.DataFrame,
    y: pd.Series,
    model_name: str = 'RandomForest',
    outer_cv: int = 5,
    inner_cv: int = 3,
    random_state: int = 42
) -> Dict[str, Any]:
    """
    Perform nested cross-validation for unbiased evaluation.

    Parameters
    ----------
    X : pd.DataFrame
        Features
    y : pd.Series
        Labels
    model_name : str
        Model to evaluate
    outer_cv : int
        Outer CV folds
    inner_cv : int
        Inner CV folds for hyperparameter tuning
    random_state : int
        Random seed

    Returns
    -------
    Dict[str, Any]
        Nested CV results
    """
    logger.info(f"Running nested CV for {model_name}")

    le = LabelEncoder()
    y_encoded = le.fit_transform(y)

    outer_scores = []
    outer_cv_split = StratifiedKFold(
        n_splits=outer_cv,
        shuffle=True,
        random_state=random_state
    )

    classifier = RNAseqClassifier(X, y, random_state=random_state)
    configs = classifier.get_model_configs()

    for fold, (train_idx, test_idx) in enumerate(outer_cv_split.split(X, y_encoded)):
        X_train_fold = X.iloc[train_idx]
        X_test_fold = X.iloc[test_idx]
        y_train_fold = y_encoded[train_idx]
        y_test_fold = y_encoded[test_idx]

        # Inner CV for hyperparameter tuning
        pipeline = classifier.create_pipeline(configs[model_name]['model'])
        inner_cv_split = StratifiedKFold(
            n_splits=inner_cv,
            shuffle=True,
            random_state=random_state
        )

        search = GridSearchCV(
            pipeline,
            configs[model_name]['params'],
            cv=inner_cv_split,
            scoring='accuracy',
            n_jobs=-1
        )
        search.fit(X_train_fold, y_train_fold)

        # Evaluate on outer test fold
        score = search.score(X_test_fold, y_test_fold)
        outer_scores.append(score)

        logger.info(f"Fold {fold + 1}: {score:.4f}")

    results = {
        'model': model_name,
        'outer_cv_scores': outer_scores,
        'mean_score': np.mean(outer_scores),
        'std_score': np.std(outer_scores),
        'ci_95': (
            np.mean(outer_scores) - 1.96 * np.std(outer_scores),
            np.mean(outer_scores) + 1.96 * np.std(outer_scores)
        )
    }

    logger.info(f"Nested CV: {results['mean_score']:.4f} +/- {results['std_score']:.4f}")

    return results


def main():
    """Test ML classifiers."""
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))

    from preprocessing.data_loader import RNAseqDataLoader
    from preprocessing.normalization import RNAseqNormalizer
    from feature_engineering.feature_selection import FeatureSelector, prepare_ml_data

    project_root = Path(__file__).parent.parent.parent
    counts_file = project_root / "data/raw/GSE313799_counts_matrix_iN.txt"
    db_path = project_root / "db/rnaseq_analysis.db"

    # Load and preprocess
    loader = RNAseqDataLoader(str(counts_file), str(db_path))
    counts = loader.load_counts()
    metadata = loader.parse_sample_metadata()

    # Filter and normalize
    keep = (counts >= 10).sum(axis=1) >= 3
    filtered = counts.loc[keep]
    normalizer = RNAseqNormalizer(filtered)
    log_cpm = normalizer.cpm(log=True)

    # Prepare data
    metadata = metadata.set_index('sample_id')
    X, y = prepare_ml_data(log_cpm, metadata, 'condition')

    # Feature selection
    selector = FeatureSelector(X, y)
    selector.variance_filter(top_n=500)
    X_selected, y_selected = selector.get_selected_data('variance')

    print(f"Selected features: {X_selected.shape}")

    # Train models
    classifier = RNAseqClassifier(X_selected, y_selected, test_size=0.2)
    results = classifier.train_all_models(tune_hyperparameters=True, cv_folds=5)

    print("\n=== Model Comparison ===")
    print(results)

    # Feature importance
    print("\n=== Top 20 Important Features (Random Forest) ===")
    importance = classifier.get_feature_importance('RandomForest')
    print(importance.head(20))

    # Save results
    results.to_csv(project_root / "results/ml_results.csv", index=False)
    importance.to_csv(project_root / "results/feature_importance.csv", index=False)

    # Save best model
    best_model = results.iloc[0]['model']
    classifier.save_model(best_model, str(project_root / f"results/{best_model}_model.joblib"))


if __name__ == "__main__":
    main()
