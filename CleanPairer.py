import pandas as pd
import os
import argparse


class CleanPairer:
    """
    Scans a variant_counts directory and generates functional selection pairs
    by matching VSVG_control (preselection) with no-antibody_control (postselection)
    files that share the same library, date, rescue batch, and replicate number.

    Expected filename format:
        {library}_{date}_{virus_batch}_{experiment_type}_{replicate}.csv
    Example:
        A_2022-07-20_rescue-2_VSVG_control_1.csv
        A_2022-07-20_rescue-2_no-antibody_control_1.csv

    Output columns:
        preselection_sample, library, virus_batch, replicate,
        postselection_sample, preselection_library_sample,
        postselection_library_sample, selection_name
    """

    def __init__(self, variant_counts_dir, output_csv):
        self.variant_counts_dir = variant_counts_dir
        self.output_csv = output_csv

    def _parse_filename(self, filename):
        name = os.path.splitext(filename)[0]
        parts = name.split("_")
        if len(parts) < 5:
            raise ValueError(f"Cannot parse filename: {filename}")
        library = parts[0]
        date = parts[1]
        virus_batch = parts[2]
        replicate = parts[-1]
        experiment_type = "_".join(parts[3:-1])
        return library, date, virus_batch, experiment_type, replicate

    def generate_pairs(self):
        files = [
            f for f in sorted(os.listdir(self.variant_counts_dir))
            if f.endswith(".csv") and not f.startswith("avg_counts")
        ]

        vsvg = {}
        no_ab = {}

        for f in files:
            try:
                library, date, virus_batch, exp_type, rep = self._parse_filename(f)
            except ValueError as e:
                print(f"Skipping: {e}")
                continue

            key = (library, date, virus_batch, rep)
            sample_with_prefix = os.path.splitext(f)[0]
            sample_no_prefix = "_".join(sample_with_prefix.split("_")[1:])

            if exp_type == "VSVG_control":
                vsvg[key] = (sample_with_prefix, sample_no_prefix)
            elif exp_type == "no-antibody_control":
                no_ab[key] = (sample_with_prefix, sample_no_prefix)

        pairs = []
        for key in sorted(vsvg):
            library, date, virus_batch, rep = key
            pre_lib, pre_no_lib = vsvg[key]
            if key in no_ab:
                post_lib, post_no_lib = no_ab[key]
                pairs.append({
                    "preselection_sample":         pre_no_lib,
                    "library":                     library,
                    "virus_batch":                 virus_batch,
                    "replicate":                   int(rep),
                    "postselection_sample":        post_no_lib,
                    "preselection_library_sample": pre_lib,
                    "postselection_library_sample": post_lib,
                    "selection_name":              f"{pre_lib}_vs_{post_no_lib}",
                })
            else:
                print(f"Warning: no matching no-antibody_control for {pre_lib}")

        return pairs

    def run(self):
        pairs = self.generate_pairs()

        if not pairs:
            raise ValueError(
                "No valid pairs found. Check that variant_counts_dir contains "
                "matching VSVG_control and no-antibody_control files."
            )

        columns = [
            "preselection_sample",
            "library",
            "virus_batch",
            "replicate",
            "postselection_sample",
            "preselection_library_sample",
            "postselection_library_sample",
            "selection_name",
        ]
        df = pd.DataFrame(pairs, columns=columns)
        print(f"Found {len(df)} valid functional pairs")

        out_dir = os.path.dirname(os.path.abspath(self.output_csv))
        os.makedirs(out_dir, exist_ok=True)

        if os.path.exists(self.output_csv):
            raise FileExistsError(
                f"Output file already exists: {self.output_csv}\n"
                "Delete it or choose a different path."
            )

        df.to_csv(self.output_csv, index=False)
        print(f"Saved: {self.output_csv}")
        return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate functional selection pairs from variant counts folder"
    )
    parser.add_argument(
        "--variant-counts-dir", required=True,
        help="Directory containing raw variant count CSVs"
    )
    parser.add_argument(
        "--output-csv", required=True,
        help="Output path for functional_selections_clean.csv"
    )
    args = parser.parse_args()

    pairer = CleanPairer(args.variant_counts_dir, args.output_csv)
    pairer.run()
