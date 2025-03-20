# After model training:
explainer = shap.DeepExplainer(model, X_train)
shap_values = explainer.shap_values(X_test)

# Save feature importance
feature_importance = pd.DataFrame({
    'Feature': features,
    'Importance': np.abs(shap_values).mean(0)
}).sort_values('Importance', ascending=False)
