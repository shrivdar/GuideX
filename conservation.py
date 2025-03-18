import numpy as np
from pathlib import Path
from typing import List, Tuple
from skbio import DNA, TabularMSA
from scipy.spatial.distance import jensenshannon
import plotly.express as px
from guidex.utils.logger import setup_logger
from Bio import AlignIO
from collections import Counter

logger = setup_logger(__name__)

class ConservationAnalyzer:
    """Enhanced conservation analysis with NaN handling and gap filtering"""
    
    def __init__(self, window_size: int = 30):
        self.window_size = window_size
        self.min_conservation = 0.7
        self.epsilon = 1e-10  # Smoothing factor for zero probabilities
        self.max_gap = 0.9    # Maximum allowed gap proportion per column

    def calculate_jsd(self, aligned_file: Path) -> Tuple[List[float], int]:
        """Calculate Jensen-Shannon Divergence with proper return structure"""
        jsd_scores = []
        valid_windows = 0
        
        try:
            alignment = AlignIO.read(aligned_file, "fasta")
            num_seqs = len(alignment)
            if num_seqs < 2:
                return [], 0  # Proper return for insufficient sequences
    
            seq_len = alignment.get_alignment_length()
            
            for i in range(0, seq_len - self.window_size + 1, self.window_size):
                window_scores = []
                window_end = i + self.window_size
                
                for j in range(i, window_end):
                    if j >= seq_len:
                        break
                    
                    col = [str(rec.seq[j]).upper() for rec in alignment]
                    if self._is_invalid_column(col):
                        continue
                    
                    # Calculate pairwise JSD for column
                    freqs = [self._safe_frequencies(nt) for nt in col]
                    column_jsd = []
                    
                    for k in range(len(freqs)):
                        for l in range(k+1, len(freqs)):
                            jsd = self._safe_jsd(freqs[k], freqs[l])
                            if not np.isnan(jsd):
                                column_jsd.append(jsd)
                    
                    if column_jsd:
                        window_scores.append(np.mean(column_jsd))
                
                if window_scores:
                    jsd_scores.append(np.mean(window_scores))
                    valid_windows += 1
                    
            return jsd_scores, valid_windows  # Correct return structure
        
        except Exception as e:
            logger.error(f"JSD calculation failed: {str(e)}")
            return [], 0  # Ensure return tuple on error

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
        # Enhanced pseudocount scheme
        counts = np.array([2.0, 2.0, 2.0, 2.0])  # Stronger pseudocounts
        nt_map = {'A':0, 'C':1, 'G':2, 'T':3}
        
        if nucleotide.upper() in nt_map:
            counts[nt_map[nucleotide.upper()]] += 3.0  # Observed nucleotide bonus
            
        # Ensure numerical stability
        total = np.sum(counts) + self.epsilon
        return np.clip(counts/total, 1e-10, 1.0)  # Prevent zero/negative values

    def _safe_jsd(self, p: np.ndarray, q: np.ndarray) -> float:
        """Numerically stable JSD calculation"""
        with np.errstate(divide='ignore', invalid='ignore'):
            p_norm = p / (np.sum(p) + self.epsilon)
            q_norm = q / (np.sum(q) + self.epsilon)
            
            if np.any(np.isnan(p_norm)) or np.any(np.isnan(q_norm)):
                return np.nan
                
            jsd = jensenshannon(p_norm, q_norm) ** 2
            return jsd if not np.isnan(jsd) else 0.0

    def _is_invalid_column(self, col: np.ndarray) -> bool:
        """Check for uninformative columns"""
        unique = len(set(col))
        return (unique < 2) or (np.mean(col == '-') > self.max_gap)

    def plot_conservation(self, scores: List[float], output_file: Path) -> None:
        """Generate interactive conservation plot with valid data structure"""
        if not scores:
            logger.warning("No scores to visualize")
            return
    
        try:
            # Create proper DataFrame structure
            positions = np.arange(len(scores))
            df = pd.DataFrame({
                'Position': positions,
                'Conservation': scores
            })
            
            # Filter NaN values
            df_clean = df.dropna()
            
            fig = px.line(
                df_clean,
                x='Position',
                y='Conservation',
                title=f"Conservation Profile (Window={self.window_size})",
                labels={'Conservation': 'JSD Score'}
            )
            
            # Add confidence band
            fig.add_trace(go.Scatter(
                x=df_clean['Position'],
                y=df_clean['Conservation'].rolling(5).mean() + 0.1,
                fill=None,
                line_color='rgba(255,0,0,0.1)',
                name='Confidence Band'
            ))
            
            output_file.parent.mkdir(parents=True, exist_ok=True)
            fig.write_html(str(output_file))
            logger.info(f"Saved conservation plot to {output_file}")
            
        except Exception as e:
            logger.error(f"Visualization failed: {str(e)}")
            raise

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
        counts = {
            'A': 0.25,  # Pseudocounts to avoid zero probabilities
            'C': 0.25,
            'G': 0.25,
            'T': 0.25
        }
        
        if nucleotide.upper() in valid_nt:
            counts[nucleotide.upper()] += 1.0
            
        total = sum(counts.values())
        return [v/total for v in counts.values()]
