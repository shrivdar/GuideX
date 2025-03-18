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
        """Guaranteed tuple return (scores, valid_window_count)"""
        try:
            alignment = AlignIO.read(aligned_file, "fasta")
            if len(alignment) < 2:
                return [], 0  # Must return tuple even when empty
                
            # Simplified working implementation
            scores = []
            valid = 0
            aln_length = alignment.get_alignment_length()
            
            for i in range(0, aln_length, self.window_size):
                window_scores = []
                for j in range(i, min(i+self.window_size, aln_length)):
                    col = [str(rec.seq[j]).upper() for rec in alignment]
                    
                    # Skip bad columns
                    if col.count('-')/len(col) > 0.8:
                        continue
                    
                    # Stable frequency calculation
                    counts = {'A':1, 'C':1, 'G':1, 'T':1}
                    for nt in col:
                        if nt in counts:
                            counts[nt] += 1
                    
                    total = sum(counts.values())
                    freqs = [v/total for v in counts.values()]
                    
                    # Pairwise JSD calculation
                    jsd_sum = 0
                    pairs = 0
                    for k in range(len(freqs)):
                        for l in range(k+1, len(freqs)):
                            jsd = jensenshannon(freqs[k], freqs[l]) ** 2
                            jsd_sum += jsd
                            pairs += 1
                    
                    if pairs > 0:
                        window_scores.append(jsd_sum/pairs)
                
                if window_scores:
                    scores.append(np.mean(window_scores))
                    valid += 1
                    
            return scores, valid  # Correct return structure
            
        except Exception as e:
            logger.error(f"Conservation analysis failed: {e}")
            return [], 0  # Always return tuple

        return jsd_scores, valid_regions

    
    def calculate_jsd(self, aligned_file: Path) -> Tuple[List[float], int]:
        """Guaranteed tuple return (scores, valid_window_count)"""
        try:
            alignment = AlignIO.read(aligned_file, "fasta")
            if len(alignment) < 2:
                return [], 0  # Proper tuple return
    
            scores = []
            valid = 0
            aln_length = alignment.get_alignment_length()
            
            for i in range(0, aln_length, self.window_size):
                window_scores = []
                for j in range(i, min(i+self.window_size, aln_length)):
                    col = [str(rec.seq[j]).upper() for rec in alignment]
                    
                    # Skip gap-heavy columns
                    if (col.count('-')/len(col)) > 0.8:
                        continue
                    
                    # Enhanced frequency calculation
                    counts = {'A': 2, 'C': 2, 'G': 2, 'T': 2}  # Stronger pseudocounts
                    for nt in col:
                        if nt in counts:
                            counts[nt] += 1
                    
                    total = sum(counts.values())
                    freqs = [v/total for v in counts.values()]
                    
                    # Calculate pairwise JSD
                    jsd_sum = 0.0
                    pairs = 0
                    for k in range(len(freqs)):
                        for l in range(k+1, len(freqs)):
                            with np.errstate(divide='ignore', invalid='ignore'):
                                jsd = jensenshannon(freqs[k], freqs[l]) ** 2
                                if not np.isnan(jsd):
                                    jsd_sum += jsd
                                    pairs += 1
                    
                    if pairs > 0:
                        window_scores.append(jsd_sum / pairs)
                
                if window_scores:
                    scores.append(np.mean(window_scores))
                    valid += 1
                    
            return scores, valid  # Correct return
    
        except Exception as e:
            logger.error(f"Conservation analysis failed: {e}")
            return [], 0  # Proper error return
    
        # REMOVE THIS LINE - CAUSES UNDEFINED VARIABLE ERRORS
        # return jsd_scores, valid_regions  

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
        """Bulletproof frequency calculation"""
        counts = np.array([2.0, 2.0, 2.0, 2.0])  # Strong pseudocounts
        nt_map = {'A':0, 'C':1, 'G':2, 'T':3}
        
        if nucleotide.upper() in nt_map:
            counts[nt_map[nucleotide.upper()]] += 3.0
            
        total = np.sum(counts)
        return np.clip(counts/total, 1e-9, 1.0)  # Prevent division issues

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
        """Simplified working plot"""
        if not scores:
            return  # Skip empty plots
        
        try:
            # Create clean dataframe
            df = pd.DataFrame({
                'Position': range(len(scores)),
                'Conservation': scores
            })
            
            # Filter NaN/Inf values
            df = df.replace([np.inf, -np.inf], np.nan).dropna()
            
            fig = px.line(
                df,
                x='Position',
                y='Conservation',
                title='Conservation Profile',
                labels={'Conservation': 'Score'}
            )
            fig.write_html(str(output_file))
            
        except Exception as e:
            logger.error(f"Plot failed: {e}")

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
