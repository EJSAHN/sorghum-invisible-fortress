import os, glob
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
# =============================================================================
# In Silico Simulation
# =============================================================================
# 1. Load Phenotypes (Spectra) & Genotype PCs from S1
# Note: Assuming S1_Phenotypes.csv contains both traits and PCs (PC1-PC10)
print("[Step] Loading real dataset from S1...")

try:
    # Load the Phenotype Sheet
    # (Update sheet_name if needed, e.g., 'S1_Phenotypes')
    df_pheno = pd.read_excel("Supplementary_Data_S1_Final_Submitted.xlsx", sheet_name="S1_Phenotypes")
    
    # Extract Spectral Traits (Target Y)
    # Assuming 'Vis_PC1', 'NIR_PC1', 'SDI' are the targets
    target_cols = ['Vis_PC1', 'NIR_PC1', 'SDI'] 
    S_data = df_pheno[target_cols].fillna(0).values
    
    # Extract Genotype PCs (Predictor X)
    # Assuming columns are named 'PC1', 'PC2', ... 'PC10'
    geno_pc_cols = [f"PC{i}" for i in range(1, 11)]
    G_data = df_pheno[geno_pc_cols].fillna(0).values
    
    print(f" -> Loaded Genotype PCs: {G_data.shape}")
    print(f" -> Loaded Spectral Targets: {S_data.shape}")

except Exception as e:
    print(f"[Error] Failed to load real data: {e}")
    # Stop execution if data is missing
    raise

# 2. Train Ridge Regression Model (Real Analysis)
print("[Step] Training Ridge Regression model (G -> S)...")
# Standardize inputs
scaler_G = StandardScaler()
scaler_S = StandardScaler()

G_scaled = scaler_G.fit_transform(G_data)
S_scaled = scaler_S.fit_transform(S_data)

# Fit Model
model = Ridge(alpha=1.0)
model.fit(G_scaled, S_scaled)
print(f" -> Model R2 Score: {model.score(G_scaled, S_scaled):.4f}")

# 3. In Silico Perturbation (Simulation)
print("[Step] Simulating evolutionary trajectory...")
# Create a copy for perturbation
G_perturbed = G_scaled.copy()

# Shift along the leading genomic axis (e.g., PC1 + 2.0 SD)
perturb_idx = 0 # Index for PC1
shift_amount = 2.0
G_perturbed[:, perturb_idx] += shift_amount

# Predict new phenotypes
S_pred_scaled = model.predict(G_perturbed)
S_pred = scaler_S.inverse_transform(S_pred_scaled)

print(" -> Simulation complete. Ready for manifold projection.")