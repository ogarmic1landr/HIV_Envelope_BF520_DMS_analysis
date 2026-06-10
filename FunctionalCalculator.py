import pandas as pd
import numpy as np
import logging


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