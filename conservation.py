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
        """Calculate Jensen-Shannon Divergence for alignment windows"""
        jsd_scores = []
        valid_windows = 0
        
        try:
            alignment = AlignIO.read(aligned_file, "fasta")
            seq_count = len(alignment)
            if seq_count < 2:
                return [], 0
                
            aln_length = alignment.get_alignment_length()
            
            for i in range(0, aln_length, self.window_size):
                window_end = min(i + self.window_size, aln_length)
                window_jsd = []
                
                # Column-wise analysis
                for col_idx in range(i, window_end):
                    column = [str(rec.seq[col_idx]).upper() for rec in alignment]
                    
                    # Skip gap-heavy columns
                    if (column.count('-') / len(column)) > self.max_gap:
                        continue
                    
                    # Calculate position frequencies
                    freqs = []
                    for nt in column:
                        counts = Counter(nt)
                        total = sum(counts.values()) + 4 * self.epsilon  # Add pseudocounts
                        freq = [
                            (counts.get('A', 0) + self.epsilon) / total,
                            (counts.get('C', 0) + self.epsilon) / total,
                            (counts.get('G', 0) + self.epsilon) / total,
                            (counts.get('T', 0) + self.epsilon) / total
                        ]
                        freqs.append(freq)
                    
                    # Calculate pairwise JSD
                    pairwise_jsd = []
                    for j in range(len(freqs)):
                        for k in range(j+1, len(freqs)):
                            jsd = jensenshannon(freqs[j], freqs[k]) ** 2
                            if not np.isnan(jsd):
                                pairwise_jsd.append(jsd)
                    
                    if pairwise_jsd:
                        window_jsd.append(np.mean(pairwise_jsd))
                
                if window_jsd:
                    jsd_scores.append(np.mean(window_jsd))
                    valid_windows += 1
                    
        except Exception as e:
            logger.error(f"JSD calculation failed: {str(e)}")
            raise
        
        return jsd_scores, valid_windows

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
        """Frequency calculation with enhanced pseudocounts"""
        counts = np.full(4, 1.0)  # Strong pseudocounts (1.0 per base)
        nt_map = {'A':0, 'C':1, 'G':2, 'T':3}
        
        if nucleotide.upper() in nt_map:
            counts[nt_map[nucleotide.upper()]] += 2.0  # Additional weight for observed base
            
        return counts / (np.sum(counts) + self.epsilon)

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
        """Generate interactive conservation plot with valid data checks"""
        if not scores:
            logger.warning("No conservation scores to visualize")
            return
    
        try:
            df = pd.DataFrame({
                "Position": np.arange(len(scores)),
                "Conservation": scores
            })
            
            fig = px.line(
                df,
                x="Position",
                y="Conservation",
                title=f"Conservation Profile (Window={self.window_size})",
                labels={"Conservation": "JSD Score"}
            )
            
            fig.update_layout(
                xaxis_rangeslider_visible=True,
                template="plotly_white"
            )
            
            output_file.parent.mkdir(exist_ok=True)
            fig.write_html(str(output_file))
            logger.info(f"Saved conservation plot to {output_file}")
            
        except Exception as e:
            logger.error(f"Visualization failed: {str(e)}")
            if 'df' in locals():
                logger.debug(f"Data summary:\n{df.describe()}")

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
