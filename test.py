import os
import pandas as pd

from FunctionalCalculator import FunctionalScoreCalculator


# =========================
# PATHS
# =========================
mapped_folder = (
    "/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/"
    "results/func_scores/merged_output_clean/"
    "mapped_output_clean"
)

selection_file = (
    "/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/"
    "results/func_scores/functional_selections_clean.csv"
)

output_folder = os.path.join(
    mapped_folder,
    "pipeline_output"
)

os.makedirs(output_folder, exist_ok=True)


# =========================
# LOAD METADATA
# =========================
selections = pd.read_csv(selection_file)

meta = selections[
    [
        "selection_name",
        "preselection_sample",
        "postselection_sample"
    ]
].rename(columns={
    "preselection_sample": "pre_sample",
    "postselection_sample": "post_sample"
})


# =========================
# INIT CALCULATOR
# =========================
calc = FunctionalScoreCalculator(
    pseudocount=0.5,
    pre_count_threshold=25
)


# =========================
# GET INPUT FILES
# =========================
files = [
    f for f in os.listdir(mapped_folder)
    if f.endswith(".csv")
    and "_pipeline" not in f
]

print(f"Found {len(files)} mapped files\n")


# =========================
# PROCESS FILES
# =========================
for file in files:

    filepath = os.path.join(
        mapped_folder,
        file
    )

    print(f"Processing: {file}")

    # -------------------------
    # Load mapped CSV
    # -------------------------
    df = pd.read_csv(filepath)

    # -------------------------
    # WT DEBUG
    # -------------------------
    wt_rows = df[
        df["n_codon_substitutions"] == 0
    ]

    print("\nWT DEBUG:")
    print(
        f"WT rows found: {len(wt_rows)}"
    )

    print(
        wt_rows[
            [
                "barcode",
                "pre_count",
                "post_count"
            ]
        ].head(10)
    )

    print(
        "\nWT TOTAL COUNTS:"
    )

    print(
        "Total WT pre_count:",
        wt_rows["pre_count"].sum()
    )

    print(
        "Total WT post_count:",
        wt_rows["post_count"].sum()
    )

    # -------------------------
    # Run calculator
    # -------------------------
    try:

        result = calc.run(df)

        # -------------------------
        # Merge metadata
        # -------------------------
        result = result.merge(
            meta,
            on="selection_name",
            how="left"
        )

        # Validate merge
        if result["pre_sample"].isna().any():

            raise ValueError(
                f"Metadata merge failed "
                f"for file: {file}"
            )

        # -------------------------
        # Rename columns
        # -------------------------
        result = result.rename(columns={
            "aa_substitutions_sequence":
            "aa_substitutions_sequential"
        })

        # -------------------------
        # Add library column
        # -------------------------
        result["library"] = (
            file.split("_")[0]
        )

        # -------------------------
        # Final column order
        # -------------------------
        final_columns = [
            "library",
            "pre_sample",
            "post_sample",
            "barcode",
            "func_score",
            "func_score_var",
            "pre_count",
            "post_count",
            "pre_count_wt",
            "post_count_wt",
            "pseudocount",
            "n_codon_substitutions",
            "aa_substitutions_sequential",
            "n_aa_substitutions",
            "aa_substitutions_reference",
            "pre_count_threshold"
        ]

        # Validate final columns
        missing_cols = [
            c for c in final_columns
            if c not in result.columns
        ]

        if missing_cols:

            raise ValueError(
                f"Missing final columns: "
                f"{missing_cols}"
            )

        result = result[
            final_columns
        ]

        # -------------------------
        # Save output
        # -------------------------
        output_path = os.path.join(
            output_folder,
            file.replace(
                ".csv",
                "_pipeline.csv"
            )
        )

        result.to_csv(
            output_path,
            index=False
        )

        print(
            f"\nSaved: {output_path}\n"
        )

    except Exception as e:

        print(
            f"\nError processing "
            f"{file}: {e}\n"
        )


print("ALL FILES PROCESSED SUCCESSFULLY")