import numpy as np
from pathlib import Path
from typing import List, Tuple
from scipy.spatial.distance import jensenshannon
from scipy.stats import entropy
from skbio import DNA, TabularMSA
from Bio import AlignIO
import pandas as pd
import plotly.express as px
from guidex.utils.logger import setup_logger

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Numerically stable conservation analysis with enhanced error handling"""
    
    def __init__(self, window_size: int = 30, max_gap: float = 0.7, pseudocount: float = 10.0):  # FIXED
        self.window_size = window_size
        self.max_gap = max_gap
        self.pseudocount = pseudocount
        self.epsilon = 1e-10

class ConservationAnalyzer:
    """Numerically stable conservation analysis with enhanced error handling"""
    
    def __init__(self, window_size: int = 30, max_gap: float = 0.7, pseudocount: float = 10.0):
        self.window_size = window_size
        self.max_gap = max_gap
        self.pseudocount = pseudocount
        self.epsilon = 1e-10

    def calculate_jsd(self, alignment_file):
        """Calculate Jensen-Shannon Divergence scores for each window"""
        align = AlignIO.read(alignment_file, "fasta")
        alignment_length = align.get_alignment_length()
        jsd_scores = []
        valid_windows = []
        
        for i in range(0, alignment_length, self.window_size):
            window_start = i
            window_end = min(i + self.window_size, alignment_length)
            window = align[:, window_start:window_end]
            window_length = window_end - window_start
            
            # Calculate gap percentage
            gap_count = sum(sum(1 for base in record.seq if base == '-') for record in window)
            total_bases = len(window) * window_length
            gap_percentage = gap_count / total_bases if total_bases > 0 else 1.0
            
            if gap_percentage > self.max_gap:
                valid_windows.append(False)
                jsd_scores.append(0.0)
                continue
                
            # Calculate frequency matrix
            freqs = np.zeros((4, window_length))
            for col in range(window_length):
                counts = {'A': self.pseudocount, 'C': self.pseudocount, 'G': self.pseudocount, 'T': self.pseudocount}
                for record in window:
                    base = record.seq[col].upper()
                    if base in counts:
                        counts[base] += 1
                    # else:  # Handle gaps and unknowns by ignoring
                
                total = sum(counts.values())
                freqs[0, col] = counts['A'] / total
                freqs[1, col] = counts['C'] / total
                freqs[2, col] = counts['G'] / total
                freqs[3, col] = counts['T'] / total
            
            # Calculate JSD
            mean_freq = np.mean(freqs, axis=1)
            jsd = 0
            for col in range(window_length):
                col_freq = freqs[:, col]
                jsd += 0.5 * (entropy(col_freq) + entropy(mean_freq) - entropy(0.5*(col_freq + mean_freq)))
            
            jsd_scores.append(jsd / window_length)
            valid_windows.append(True)
        return np.array(jsd_scores), valid_windows

    def _safe_jsd(self, p: List[float], q: List[float]) -> float:
        """Numerically stable JSD implementation"""
        p = np.asarray(p)
        q = np.asarray(q)
        
        # Add epsilon to prevent division by zero
        p = np.clip(p, 1e-12, 1.0)
        q = np.clip(q, 1e-12, 1.0)
        
        m = 0.5 * (p + q)
        return 0.5 * (entropy(p, m) + entropy(q, m))

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Robust plotting with NaN handling"""
        try:
            if not scores:
                logger.warning("No scores to plot")
                return

            # Create clean dataframe with position index
            df = pd.DataFrame({
                'Position': np.arange(len(scores)),
                'Conservation': np.nan_to_num(scores, nan=0.0)
            })
            
            # Generate interactive plot
            fig = px.line(
                df,
                x='Position',
                y='Conservation',
                title=f'Conservation Profile (Window Size: {self.window_size}bp)',
                labels={'Conservation': 'Jensen-Shannon Distance'},
                template='plotly_white'
            )
            
            fig.write_html(str(output_file))
            logger.info(f"Saved conservation plot to {output_file}")

        except Exception as e:
            logger.error(f"Visualization error: {str(e)}")
            if hasattr(e, 'message'):
                logger.debug(f"Error details: {e.message}")
    
    def _load_and_filter_alignment(self, path: Path) -> TabularMSA:
        """Load MSA and filter gap-heavy columns"""
        msa = TabularMSA.read(str(path), constructor=DNA)
        
        # Convert to gap-filtered numpy array
        seq_array = np.array([list(str(seq)) for seq in msa])
        gap_freq = np.mean(seq_array == '-', axis=0)
        valid_cols = gap_freq <= self.max_gap
        filtered = seq_array[:, valid_cols]
        
        logger.info(f"Filtered {len(gap_freq)-sum(valid_cols)}/{len(gap_freq)} gap-heavy columns")
        return TabularMSA([DNA(''.join(row)) for row in filtered])

    def _safe_windowed_analysis(self, msa: TabularMSA) -> Tuple[List[float], int]:
        """JSD calculation with numerical stability checks"""
        seq_array = np.array([[nt for nt in str(seq)] for seq in msa])
        num_seqs, seq_len = seq_array.shape
        scores = []
        valid_windows = 0

        for i in range(0, seq_len - self.window_size + 1):
            window_scores = []
            for j in range(i, i + self.window_size):
                if j >= seq_len:
                    break
                
                col = seq_array[:, j]
                if self._is_invalid_column(col):
                    continue

                # Calculate pairwise JSD with smoothing
                jsds = []
                for k in range(num_seqs):
                    for l in range(k+1, num_seqs):
                        p = self._safe_frequencies(col[k])
                        q = self._safe_frequencies(col[l])
                        jsd = self._safe_jsd(p, q)
                        if not np.isnan(jsd):
                            jsds.append(jsd)
                
                if jsds:
                    window_scores.append(np.mean(jsds))
                    valid_windows += 1
            
            if window_scores:
                scores.extend(window_scores)

        return scores, valid_windows

    def _safe_frequencies(self, nucleotide: str) -> np.ndarray:
        """Robust frequency calculation with enhanced pseudocounts"""
        # Changed from fixed 5.0 to use instance's pseudocount
        counts = np.array([self.pseudocount] * 4)  # A, C, G, T
        nt_map = {'A':0, 'C':1, 'G':2, 'T':3}
        
        if nucleotide.upper() in nt_map:
            counts[nt_map[nucleotide.upper()]] += 2.0  # Match calculate_jsd's increment
            
        total = np.sum(counts) + 1e-12
        return np.clip(counts/total, 1e-12, 1.0)

    def _is_invalid_column(self, col: np.ndarray) -> bool:
        """Check for uninformative columns"""
        # Changed from set() to np.unique for numpy arrays
        return (len(np.unique(col)) < 2) or (np.mean(col == '-') > self.max_gap)

    def _load_alignment(self, path: Path) -> TabularMSA:
        """Load and validate alignment file"""
        logger.debug(f"Loading alignment from {path}")
        msa = TabularMSA.read(str(path), constructor=DNA)
        
        # Validate alignment content
        if len(msa) == 0:
            raise ValueError("Empty alignment file")
        if any(len(seq) != len(msa[0]) for seq in msa):
            raise ValueError("Inconsistent sequence lengths in alignment")
            
        return msa

    def _windowed_analysis(self, msa: TabularMSA) -> List[float]:
        """Sliding window conservation scoring"""
        scores = []
        seq_length = len(msa[0])
        
        for i in range(0, seq_length - self.window_size + 1):
            window_scores = []
            for j in range(i, min(i+self.window_size, seq_length)):
                position_scores = []
                for seq1 in msa:
                    for seq2 in msa:
                        if seq1 != seq2:
                            p = self._position_frequencies(seq1[j])
                            q = self._position_frequencies(seq2[j])
                            position_scores.append(jensenshannon(p, q) ** 2)
                window_scores.append(np.mean(position_scores) if position_scores else 0.0)
            scores.extend(window_scores)
            
        return scores

    def _position_frequencies(self, nucleotide: str) -> List[float]:
        """Calculate normalized nucleotide frequencies at a position"""
        valid_nt = {'A', 'C', 'G', 'T'}
        # Use instance's pseudocount instead of fixed 0.25
        counts = {
            'A': self.pseudocount/40,  # Scaled for compatibility
            'C': self.pseudocount/40,
            'G': self.pseudocount/40,
            'T': self.pseudocount/40
        }
        
        if nucleotide.upper() in valid_nt:
            counts[nucleotide.upper()] += 1.0
            
        total = sum(counts.values())
        return [v/total for v in counts.values()]
