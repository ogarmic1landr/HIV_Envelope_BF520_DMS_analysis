import pandas as pd
import numpy as np
import os
import re
import logging




class MutationMapper:
    def __init__(self, input_folder, map_file, output_folder="mapped_output", suffix="_mapped"):
        self.input_folder = input_folder
        self.map_file = map_file
        self.suffix = suffix
        self.output_folder = os.path.join(input_folder, output_folder)
        os.makedirs(self.output_folder, exist_ok=True)

        # Set up a logger for this instance
        self.logger = logging.getLogger(f"{__name__}.{os.path.basename(self.output_folder)}")
        # Avoid adding handlers repeatedly if called multiple times
        if not self.logger.handlers:
            log_path = os.path.join(self.output_folder, "mutation_mapping.log")
            fh = logging.FileHandler(log_path)
            fh.setLevel(logging.INFO)
            ch = logging.StreamHandler()
            ch.setLevel(logging.INFO)
            formatter = logging.Formatter("%(levelname)s: %(message)s")
            fh.setFormatter(formatter)
            ch.setFormatter(formatter)
            self.logger.addHandler(fh)
            self.logger.addHandler(ch)
            self.logger.setLevel(logging.INFO)

        self.logger.info("Initializing MutationMapper...")
        self.map_df = pd.read_csv(self.map_file)

        # Validate FIRST (before using columns)
        required_map_cols = ["sequential_site", "reference_site"]
        missing = [c for c in required_map_cols if c not in self.map_df.columns]
        if missing:
            raise ValueError(f"Mapping file missing columns: {missing}")

        #  Then enforce types
        self.map_df["sequential_site"] = self.map_df["sequential_site"].astype(int)
        self.map_df["reference_site"] = self.map_df["reference_site"].astype(str)

        # Then build mapping
        self.mapping_dict = dict(
            zip(self.map_df["sequential_site"], self.map_df["reference_site"])
        )
        

    
    
    def count_substitutions(self, mutation_string):
        if pd.isna(mutation_string) or str(mutation_string).strip() == "":
            return 0
        s = re.sub(r"\s*-\s*", " ", str(mutation_string)).upper()
        return len([m for m in s.split() if m != "-"])
    

    
    
    def map_aa_substitutions(self, aa_string):
        if pd.isna(aa_string) or str(aa_string).strip() == "":
            return ""
        # Replace dashes with spaces and uppercase for consistency
        aa_string = str(aa_string).upper()
        mapped = []
        for mut in aa_string.split():
            if mut == "-":
                continue
            match = re.match(r"^([A-Z\*])(\d+)([A-Z\*\-]?)$", mut)
            if match:
                wt, pos, mut_aa = match.groups()
                pos = int(pos)
                new_pos = self.mapping_dict.get(pos, pos)
                if pos != new_pos:
                    self.logger.debug(f"Mapped {pos} -> {new_pos}")
                mapped.append(f"{wt}{new_pos}{mut_aa}")
            else:
                self.logger.warning(f"Malformed mutation: {mut}")
                mapped.append(mut)
        return " ".join(mapped)
    

    
    
    def process_file(self, filepath):
        self.logger.info(f"Processing file: {filepath}")
        df = pd.read_csv(filepath)
        required_cols = [
            "selection_name", "barcode", "pre_count", "post_count",
            "codon_substitutions", "aa_substitutions",
            "pre_variant_call_support", "post_variant_call_support"
        ]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(f"{filepath} missing columns: {missing}")

        df["aa_substitutions_sequence"] = df["aa_substitutions"]
        df["n_aa_substitutions"] = df["aa_substitutions_sequence"].apply(self.count_substitutions)
        df["n_codon_substitutions"] = df["codon_substitutions"].apply(self.count_substitutions)
        df["aa_substitutions_reference"] = df["aa_substitutions_sequence"].apply(self.map_aa_substitutions)
        df = df.drop(columns=["aa_substitutions"])

        base = os.path.splitext(os.path.basename(filepath))[0]
        output_filename = f"{base}{self.suffix}.csv"
        output_path = os.path.join(self.output_folder, output_filename)

        # Optional: prevent overwriting
        if os.path.exists(output_path):
            self.logger.warning(f"Output file exists, overwriting: {output_filename}")
        df.to_csv(output_path, index=False)
        self.logger.info(f"Saved: {output_path}")



    def run_all(self):
        files = [f for f in os.listdir(self.input_folder) if f.endswith(".csv")]
        if not files:
            self.logger.warning("No CSV files found in input folder")
            return
        self.logger.info(f"Found {len(files)} files to process")
        for file in files:
            self.process_file(os.path.join(self.input_folder, file))
        self.logger.info("All files processed successfully")




if __name__ == "__main__":
    # Set up root logger once
    #logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    mapper = MutationMapper(
        input_folder="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/func_scores/merged_output_clean",
        map_file="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/site_numbering/site_numbering_map.csv",
        output_folder="mapped_output_clean",
        suffix="_mapped"
    )
    mapper.run_all()