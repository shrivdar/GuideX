import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import OneHotEncoder
from .grna_designer import Cas13gRNADesigner  # Inherit base constraints

class Cas13Optimizer(nn.Module):
    """Neural optimizer combining DeepCas13 architecture with rule-based constraints"""
    
    def __init__(self, designer: Cas13gRNADesigner):
        super().__init__()
        self.designer = designer
        self.embed_size = 64
        
        # Neural network components
        self.conv = nn.Sequential(
            nn.Conv1d(4, 32, kernel_size=5),
            nn.ReLU(),
            nn.MaxPool1d(2),
            nn.Dropout(0.2)
        )
        
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

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        # Input shape: (batch_size, seq_len, 4)
        x = x.permute(0, 2, 1)  # Conv1d expects channels first
        x = self.conv(x)
        x = x.permute(0, 2, 1)
        
        # LSTM processing
        lstm_out, _ = self.lstm(x)
        
        # Attention mechanism
        attn_out, _ = self.attention(
            lstm_out, lstm_out, lstm_out
        )
        
        # Combine features
        combined = torch.cat([lstm_out, attn_out], dim=-1)
        pooled = torch.mean(combined, dim=1)
        
        return self.classifier(pooled)

    def optimize(self, spacer: str, n_iter: int = 50) -> str:
        """Optimize spacer sequence using gradient ascent"""
        seq = spacer.upper()
        best_score = -np.inf
        best_seq = seq
        
        for _ in range(n_iter):
            # Convert to one-hot encoding
            seq_tensor = self._sequence_to_tensor(seq)
            seq_tensor.requires_grad_(True)
            
            # Predict score
            score = self(seq_tensor.unsqueeze(0))
            
            # Backpropagate
            self.zero_grad()
            score.backward()
            
            # Get gradient-guided mutation
            with torch.no_grad():
                grad = seq_tensor.grad.numpy()
                mutated = self._apply_gradient_mutation(seq, grad)
                
                # Validate and score
                if self._is_valid(mutated):
                    new_score = self(self._sequence_to_tensor(mutated).unsqueeze(0))
                    if new_score > best_score:
                        best_score = new_score
                        best_seq = mutated
                        
        return best_seq if best_score > self(best_seq) else spacer

    def _is_valid(self, seq: str) -> bool:
        """Check against designer constraints"""
        return self.designer._validate_basic(
            seq, 
            self.designer._calculate_gc(seq)
        )

    def _sequence_to_tensor(self, seq: str) -> torch.Tensor:
        """One-hot encode DNA sequence"""
        encoder = OneHotEncoder(categories=[['A','C','G','T']], sparse=False)
        encoded = encoder.fit_transform(list(seq))
        return torch.tensor(encoded.T, dtype=torch.float32)

    def _apply_gradient_mutation(self, seq: str, grad: np.ndarray) -> str:
        """Apply gradient-guided nucleotide substitutions"""
        new_seq = []
        for i, nt in enumerate(seq):
            nt_grad = grad[:, i]
            candidates = []
            
            for mut in ['A','C','G','T']:
                if mut == nt:
                    continue
                # Calculate gradient similarity
                mut_vec = self._sequence_to_tensor(mut).numpy()
                score = np.dot(nt_grad, mut_vec)
                candidates.append((score, mut))
                
            # Select top beneficial mutation
            if candidates:
                _, best_mut = max(candidates)
                new_seq.append(best_mut)
            else:
                new_seq.append(nt)
                
        return ''.join(new_seq)

class OptimizationDataset(Dataset):
    """Adapter for loading gRNA activity datasets"""
    def __init__(self, sequences: List[str], activities: List[float]):
        self.encoder = OneHotEncoder(categories=[['A','C','G','T']], sparse=False)
        self.X = [self._encode(s) for s in sequences]
        self.y = activities

    def _encode(self, seq: str) -> torch.Tensor:
        encoded = self.encoder.fit_transform(list(seq))
        return torch.tensor(encoded.T, dtype=torch.float32)

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        return self.X[idx], self.y[idx]

def train_optimizer(model: Cas13Optimizer, 
                   train_data: DataLoader,
                   epochs: int = 100,
                   lr: float = 0.001):
    """Training loop with constraint-aware loss"""
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr)
    loss_fn = nn.MSELoss()
    
    for epoch in range(epochs):
        for X, y in train_data:
            pred = model(X)
            loss = loss_fn(pred, y)
            
            # Add constraint regularization
            validity_loss = torch.mean(
                torch.tensor([model._is_valid(seq) for seq in X], 
                dtype=torch.float32)
            )
            loss += 0.1 * (1 - validity_loss)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
