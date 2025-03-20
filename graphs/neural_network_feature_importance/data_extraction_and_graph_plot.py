import numpy as np
import matplotlib.pyplot as plt
import shap
import tensorflow as tf
from typing import Dict, List

class NNFeatureExplainer:
    def __init__(self, model_path: str, background_data: np.ndarray):
        """
        Initialize SHAP explainer with trained model
        Args:
            model_path: Path to saved Keras model
            background_data: Representative sample for SHAP background (50-100 samples)
        """
        self.model = tf.keras.models.load_model(model_path)
        self.explainer = shap.DeepExplainer(self.model, background_data)
        
    def calculate_shap_values(self, X: np.ndarray) -> np.ndarray:
        """Calculate SHAP values for input samples"""
        return self.explainer.shap_values(X)

def visualize_feature_importance(shap_values: np.ndarray,
                               feature_names: List[str],
                               sample_index: int = 0,
                               top_n: int = 15) -> plt.Figure:
    """
    Create horizontal bar plot of SHAP feature importance
    Args:
        shap_values: Matrix of SHAP values from NNFeatureExplainer
        feature_names: List of feature names
        sample_index: Which sample to visualize (0 for global importance)
        top_n: Number of top features to display
    """
    with plt.rc_context(POSTER_STYLE):
        # Aggregate global importance if sample_index is None
        if sample_index == 'global':
            shap_global = np.abs(shap_values).mean(0)
            values = shap_global.mean(0)
            title = "Global Feature Importance"
        else:
            values = shap_values[sample_index]
            title = f"Feature Importance (Sample {sample_index})"

        # Sort features
        indices = np.argsort(np.abs(values))[::-1][:top_n]
        sorted_names = [feature_names[i] for i in indices]
        sorted_values = values[indices]
        colors = ['#2c7bb6' if v > 0 else '#d7191c' for v in sorted_values]

        fig, ax = plt.subplots(figsize=(12, 8))
        bars = ax.barh(sorted_names, sorted_values, color=colors, ec='k', height=0.8)

        # Value annotations
        for bar in bars:
            width = bar.get_width()
            ax.text(width/2, bar.get_y() + bar.get_height()/2,
                   f"{width:.2f}", ha='center', va='center',
                   color='white', fontweight='bold', fontsize=12)

        # Formatting
        ax.set_xlabel("SHAP Value Impact", labelpad=10)
        ax.set_title(f"Neural Network {title}", pad=20)
        ax.grid(axis='x', alpha=0.3)
        
        # Panel label
        ax.text(0.01, 0.97, 'D', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')

        # Annotation box
        ax.text(0.95, 0.15, 
               "SHAP values quantify feature contribution\n"
               "to model predictions:\n"
               "• Blue = Positive impact\n"
               "• Red = Negative impact",
               transform=ax.transAxes, ha='right', va='bottom',
               bbox=dict(facecolor='white', alpha=0.8))

        plt.tight_layout()
        return fig

# Usage in analysis pipeline
def nn_feature_analysis_pipeline():
    # Load data and model
    X_train = np.load("data/processed/X_train.npy")  # Preprocessed training data
    feature_names = load_feature_names("data/processed/features.json")  # Your feature names
    
    # Initialize explainer with background data
    explainer = NNFeatureExplainer(
        model_path="models/crispr_nn_model.h5",
        background_data=X_train[:100]
    )
    
    # Calculate SHAP values for first 100 samples
    shap_values = explainer.calculate_shap_values(X_train[:100])
    
    # Generate global importance plot
    fig_global = visualize_feature_importance(
        shap_values=shap_values,
        feature_names=feature_names,
        sample_index='global'
    )
    fig_global.savefig("results/nn_feature_importance_global.png")
    
    # Generate example individual explanation
    fig_sample = visualize_feature_importance(
        shap_values=shap_values,
        feature_names=feature_names,
        sample_index=42
    )
    fig_sample.savefig("results/nn_feature_importance_sample42.png")
def plot_feature_importance(importance: Dict[str, float], 
                          top_n: int = 15):
    with plt.rc_context(POSTER_STYLE):
        features = sorted(importance.items(), key=lambda x: x[1], reverse=True)[:top_n]
        names = [f[0] for f in features]
        values = [f[1] for f in features]
        
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Colormap based on sign
        colors = ['#2c7bb6' if v >=0 else '#d7191c' for v in values]
        
        bars = ax.barh(names, values, color=colors, ec='k', height=0.8)
        
        # Value annotations
        for bar in bars:
            width = bar.get_width()
            ax.text(width/2, bar.get_y() + bar.get_height()/2,
                   f"{width:.2f}", ha='center', va='center',
                   color='white', fontweight='bold', fontsize=12)
            
        # Formatting
        ax.set_xlabel("Feature Impact Score", labelpad=10)
        ax.set_title("CRISPR Efficacy Predictor Feature Importance", pad=20)
        ax.grid(axis='x', alpha=0.3)
        
        # Panel label
        ax.text(0.01, 0.97, 'D', transform=ax.transAxes,
               fontsize=24, weight='bold', va='top')
        
        # Annotation box
        annotation_text = (
            "Features derived from:\n"
            "- Sequence composition\n"
            "- Thermodynamic modeling\n"
            "- Conservation analysis\n"
            "- Deep learning model"
        )
        ax.text(0.95, 0.25, annotation_text, transform=ax.transAxes,
               ha='right', va='bottom', fontsize=12,
               bbox=dict(facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        return fig
