from typing import List
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import OneHotEncoder
from guidex.core.grna_designer import Cas13gRNADesigner
from typing import Tuple, Dict, Optional

class Cas13Optimizer(nn.Module):
    """Neural optimizer combining DeepCas13 architecture with rule-based constraints"""
    
    def __init__(self, designer: Cas13gRNADesigner):
        super().__init__()
        self.designer = designer
        self.embed_size = 64
        
        # Initialize encoder once for all sequences
        self.encoder = OneHotEncoder(
            categories=[['A','C','G','T']], 
            sparse_output=False  # Fixed parameter name
        )
        self.encoder.fit([['A'], ['C'], ['G'], ['T']])  # Fit once during init
        
        # Neural network components
        self.conv = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=5),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(0.2)
        )
        
        # Rest of network remains the same...
        self.lstm = nn.LSTM(
            input_size=32, 
            hidden_size=self.embed_size,
            bidirectional=True
        )
        
        self.attention = nn.MultiheadAttention(
            embed_dim=self.embed_size*2,
            num_heads=4
        )
        
        self.classifier = nn.Sequential(
            nn.Linear(self.embed_size*4, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
            nn.Sigmoid()
        )


    def _sequence_to_tensor(self, seq: str) -> torch.Tensor:
        """Properly encode DNA sequence of length N to (4, N) tensor"""
        # Reshape sequence to (n_samples, n_features) = (28, 1)
        nucleotides = np.array(list(seq)).reshape(-1, 1)
        
        # Transform using pre-fit encoder
        encoded = self.encoder.transform(nucleotides)
        
        # Convert to (channels, length) format
        return torch.tensor(encoded.T, dtype=torch.float32)

    def optimize(self, spacer: str, n_iter: int = 50) -> str:
        """Optimize spacer sequence using gradient ascent"""
        seq = spacer.upper()
        best_score = -np.inf
        best_seq = seq
        
        for _ in range(n_iter):
            seq_tensor = self._sequence_to_tensor(seq)
            seq_tensor.requires_grad_(True)
            
            score = self(seq_tensor.unsqueeze(0))
            
            self.zero_grad()
            score.backward()
            
            with torch.no_grad():
                grad = seq_tensor.grad.numpy()
                mutated = self._apply_gradient_mutation(seq, grad)
                
                if self._is_valid(mutated):
                    # Compare scores properly
                    mutated_tensor = self._sequence_to_tensor(mutated).unsqueeze(0)
                    new_score = self(mutated_tensor)
                    if new_score > best_score:
                        best_score = new_score
                        best_seq = mutated
                        
        # Final validation with proper tensor conversion
        best_tensor = self._sequence_to_tensor(best_seq).unsqueeze(0)
        return best_seq if best_score > self(best_tensor) else spacer

class OptimizationDataset(Dataset):
    """Adapter for loading gRNA activity datasets"""
    def __init__(self, sequences: list[str], activities: list[float]):
        self.encoder = OneHotEncoder(
            categories=[['A','C','G','T']], 
            sparse_output=False  # Fixed parameter
        )
        self.encoder.fit([['A'], ['C'], ['G'], ['T']])
        
        self.X = [self._encode(s) for s in sequences]
        self.y = activities
        self.sequences = sequences  # Store originals for validation

    def __getitem__(self, idx: int) -> Tuple[torch.Tensor, float, str]:
        """Return (encoded, activity, original_seq)"""
        return self.X[idx], self.y[idx], self.sequences[idx]

def train_optimizer(model: Cas13Optimizer, 
                   train_data: DataLoader,
                   epochs: int = 100,
                   lr: float = 0.001):
    """Training loop with proper sequence validation"""
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()
    
    for epoch in range(epochs):
        for X, y, seqs in train_data:  # Now gets original sequences
            pred = model(X)
            loss = loss_fn(pred, y)
            
            # Validate using original sequences
            validity_loss = torch.mean(
                torch.tensor(
                    [float(model._is_valid(seq)) for seq in seqs],
                    dtype=torch.float32
                )
            )
            loss += 0.1 * (1 - validity_loss)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
