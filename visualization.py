from pathlib import Path
import logging

logger = logging.getLogger(__name__)

class PlotTracker:
    @staticmethod
    def verify_plot(grna_sequence: str, results_dir: Path = Path("results")) -> dict:
        """Check plot existence and location"""
        spacer_hash = abs(hash(grna_sequence)) % 1000
        plot_dir = results_dir / "off_targets" / f"{grna_sequence[:8]}_{spacer_hash:03d}"
        plot_path = plot_dir / "mismatch_distribution.png"
        
        return {
            "exists": plot_path.exists(),
            "path": str(plot_path),
            "directory": str(plot_dir)
        }

    @classmethod
    def report_all_plots(cls, grnas: list, results_dir: Path):
        """Generate full plot report"""
        return {grna.sequence: cls.verify_plot(grna.sequence, results_dir) 
                for grna in grnas}
