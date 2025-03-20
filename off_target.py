nano 













    def analyze(self, spacer: str) -> List[OffTarget]:
        """Main analysis workflow"""
        spacer = self._normalize_spacer(spacer)
        if not self._validate_spacer(spacer):
            return []
            
        spacer_dir = self._create_output_dir(spacer)
        targets = []
        
        # Search both forward and reverse strands
        for strand, seq in [('+', spacer), 
                          ('-', self._reverse_complement(spacer))]:
            for chrom, chrom_seq in self.genome.items():
                targets += self._find_matches(
                    search_seq=seq,
                    chromosome=chrom,
                    sequence=chrom_seq,
                    strand=strand
                )
        
        if targets:
            self._save_results(spacer, targets, spacer_dir)
            self._generate_plots(spacer, targets, spacer_dir)
            
        self._save_metadata(spacer, targets, spacer_dir)
        return targets

    def _save_results(self, spacer: str, targets: List[OffTarget], 
                    output_dir: Path):
        """Save results to CSV and JSON"""
        # CSV output
        df = pd.DataFrame([self._target_to_dict(t) for t in targets])
        df.to_csv(output_dir / "offtargets.csv", index=False)
        
        # JSON summary
        summary = {
            "spacer": spacer,
            "total_offtargets": len(targets),
            "mismatch_distribution": self._get_mismatch_counts(targets),
            "top_hits": [self._target_to_dict(t) for t in targets[:5]]
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2))

    def _save_metadata(self, spacer: str, targets: List[OffTarget], 
                     output_dir: Path):
        """Save analysis metadata"""
        metadata = {
            "spacer": spacer,
            "valid": self._validate_spacer(spacer),
            "offtarget_count": len(targets),
            "plot_exists": (output_dir / "mismatch_distribution.png").exists(),
            "analysis_dir": str(output_dir)
        }
        (output_dir / "metadata.json").write_text(json.dumps(metadata))
