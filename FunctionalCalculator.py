import pandas as pd
import numpy as np
import logging
import glob
import os
import argparse


class FunctionalScoreCalculator:

    def __init__(
        self,
        pseudocount=0.5,
        min_preselection_counts=20,
        min_preselection_frac=1e-6,
        pre_count_col="pre_count",
        post_count_col="post_count",
        codon_subs_col="n_codon_substitutions",
        selection_col="selection_name"
    ):

        self.pseudocount = pseudocount
        self.min_preselection_counts = (min_preselection_counts)

        self.min_preselection_frac = (min_preselection_frac)

        # Column configuration
        self.pre_col = pre_count_col
        self.post_col = post_count_col
        self.codon_subs_col = codon_subs_col
        self.sel_col = selection_col

        # Logger
        self.logger = logging.getLogger(
            self.__class__.__name__
        )

        if not self.logger.handlers:

            handler = logging.StreamHandler()

            formatter = logging.Formatter(
                "%(levelname)s: %(message)s"
            )

            handler.setFormatter(formatter)

            self.logger.addHandler(handler)

            self.logger.setLevel(logging.INFO)

    # -----------------------------------
    # Compute WT from raw data
    # -----------------------------------
    def compute_wt(self, df):

        self.logger.info(
            "Computing WT counts..."
        )

        # TRUE WT:
        # exact codon match to reference
        wt_mask = (
            df[self.codon_subs_col] == 0
        )

        wt_df = df[wt_mask]

        if wt_df.empty:

            raise ValueError(
                f"No WT rows found "
                f"({self.codon_subs_col} == 0)"
            )

        self.logger.info(
            f"WT rows found: {len(wt_df)}"
        )

        # Sum WT counts per selection
        wt_summary = (
            wt_df
            .groupby(self.sel_col)
            .agg({
                self.pre_col: "sum",
                self.post_col: "sum"
            })
            .rename(columns={
                self.pre_col: "pre_count_wt",
                self.post_col: "post_count_wt"
            })
            .reset_index()
        )

        self.logger.info(
            f"WT computed for "
            f"{len(wt_summary)} selections"
        )

        return wt_summary

    # -----------------------------------
    # Compute functional score
    # -----------------------------------
    def compute_scores(self, df):

        self.logger.info(
            "Computing functional scores..."
        )

        p = self.pseudocount

        # Functional score
        df["func_score"] = np.log2(
            (
                (df[self.post_col] + p) /
                (df["post_count_wt"] + p)
            ) /
            (
                (df[self.pre_col] + p) /
                (df["pre_count_wt"] + p)
            )
        )

        # Functional score variance
        df["func_score_var"] = (
            (1 / (np.log(2) ** 2)) *
            (
                1 / (df[self.post_col] + p) +
                1 / (df[self.pre_col] + p) +
                1 / (df["post_count_wt"] + p) +
                1 / (df["pre_count_wt"] + p)
            )
        )

        return df

    # -----------------------------------
    # Main pipeline
    # -----------------------------------
    def run(self, df):

        self.logger.info(
            "Starting Functional Score Calculation..."
        )

        df = df.copy()

        # Validate required columns
        required_cols = [
            self.sel_col,
            self.pre_col,
            self.post_col,
            self.codon_subs_col
        ]

        missing = [
            c for c in required_cols
            if c not in df.columns
        ]

        if missing:

            raise ValueError(
                f"Missing required columns: {missing}"
            )

        # -----------------------------------
        # Step 1: Compute WT
        # -----------------------------------
        wt_summary = self.compute_wt(df)

        # -----------------------------------
        # Step 2: Merge WT
        # -----------------------------------
        df = df.merge(
            wt_summary,
            on=self.sel_col,
            how="left"
        )

        # -----------------------------------
        # Step 3: Missing WT
        # -----------------------------------
        missing_mask = (
            df["pre_count_wt"].isna() |
            df["post_count_wt"].isna()
        )

        if missing_mask.any():

            missing = (
                df.loc[
                    missing_mask,
                    self.sel_col
                ]
                .unique()
            )

            raise ValueError(
                f"WT computation failed.\n"
                f"Missing WT values for "
                f"selection(s): "
                f"{list(sorted(missing))}"
            )

        # -----------------------------------
        # Step 4: Invalid WT
        # -----------------------------------
        invalid_mask = (
            (df["pre_count_wt"] <= 0) |
            (df["post_count_wt"] <= 0)
        )

        if invalid_mask.any():

            bad = (
                df.loc[
                    invalid_mask,
                    self.sel_col
                ]
                .unique()
            )

            raise ValueError(
                f"Invalid WT values detected.\n"
                f"WT counts are zero or negative "
                f"for selection(s): "
                f"{list(sorted(bad))}"
            )

        total_preselection_counts = (
            df[self.pre_col].sum()
        )


        dynamic_threshold = round(
            total_preselection_counts * 
            self.min_preselection_frac
        )


        dynamic_threshold = max(
            self.min_preselection_counts,
            dynamic_threshold
        )


        self.logger.info(
            f"Dynamic pre-count threshold: "
            f"{dynamic_threshold}"
        )



        # -----------------------------------
        # Step 5: Compute scores
        # -----------------------------------
        df = self.compute_scores(df)

        # -----------------------------------
        # Step 6: Metadata
        # -----------------------------------
        df["pseudocount"] = (
            self.pseudocount
        )

        df["pre_count_threshold"] = (
            dynamic_threshold
        )



        self.logger.info(
            "Functional score calculation complete"
        )

        return df

    # -----------------------------------
    # Process all mapped files in a folder
    # -----------------------------------
    def run_all(self, input_dir, output_dir, selections_file=None):

        os.makedirs(output_dir, exist_ok=True)

        # Build pre_sample / post_sample lookup from selections CSV
        meta = None
        if selections_file:
            sel = pd.read_csv(selections_file)
            meta = sel[
                ["selection_name", "preselection_sample", "postselection_sample"]
            ].rename(columns={
                "preselection_sample": "pre_sample",
                "postselection_sample": "post_sample",
            })

        files = sorted(glob.glob(os.path.join(input_dir, "*.csv")))

        if not files:
            raise FileNotFoundError(
                f"No CSV files found in {input_dir}"
            )

        self.logger.info(f"Found {len(files)} files to process")

        final_columns = [
            "library", "pre_sample", "post_sample", "barcode",
            "func_score", "func_score_var",
            "pre_count", "post_count", "pre_count_wt", "post_count_wt",
            "pseudocount", "n_codon_substitutions",
            "aa_substitutions_sequential", "n_aa_substitutions",
            "aa_substitutions_reference", "pre_count_threshold",
        ]

        for filepath in files:

            filename = os.path.basename(filepath)
            self.logger.info(f"Processing: {filename}")

            df = pd.read_csv(filepath)
            result = self.run(df)

            result["library"] = filename.split("_")[0]

            result = result.rename(columns={
                "aa_substitutions_sequence": "aa_substitutions_sequential"
            })

            if meta is not None:
                result = result.merge(meta, on="selection_name", how="left")

            missing = [c for c in final_columns if c not in result.columns]
            if missing:
                raise ValueError(f"Missing final columns in {filename}: {missing}")

            result = result[final_columns]

            output_name = filename.replace(".csv", "_func_score.csv")
            output_path = os.path.join(output_dir, output_name)
            result.to_csv(output_path, index=False)
            self.logger.info(f"Saved: {output_path}")

        self.logger.info("All files processed successfully")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute functional scores from mapped variant count CSVs"
    )
    parser.add_argument(
        "--input-dir", required=True,
        help="Directory containing mapped CSV files (output of MutationMapper)"
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Output directory for *_func_score.csv files"
    )
    parser.add_argument(
        "--selections", default=None,
        help="Path to functional_selections_clean.csv for pre_sample/post_sample metadata"
    )
    parser.add_argument(
        "--pseudocount", type=float, default=0.5
    )
    parser.add_argument(
        "--min-preselection-counts", type=int, default=20
    )
    parser.add_argument(
        "--min-preselection-frac", type=float, default=1e-6
    )
    args = parser.parse_args()

    calc = FunctionalScoreCalculator(
        pseudocount=args.pseudocount,
        min_preselection_counts=args.min_preselection_counts,
        min_preselection_frac=args.min_preselection_frac,
    )
    calc.run_all(args.input_dir, args.output_dir, selections_file=args.selections)