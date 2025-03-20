import pandas as pd

# Mock CSV structure:
# gRNA_ID,Predicted_Efficiency,Validated_Efficiency,Selected
# gRNA-RdRp1,87.2,85.4,True
# gRNA-001,72.1,68.3,False

df = pd.read_csv("results/grna_scores.csv")
selected_grnas = df[df['Selected'] == True]
