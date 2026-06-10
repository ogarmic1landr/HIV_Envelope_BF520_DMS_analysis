import pandas as pd
import os


class CleanPairer:
    def __init__(self, input_csv, output_csv):
        self.input_csv = input_csv
        self.output_csv = output_csv

    def parse_sample(self, sample):
        """
        Extract prefix, date, and replicate from sample name
        Example:
        A_2022-09-01_rescue-3_VSVG_control_1
        """
        parts = sample.split("_")

        prefix = parts[0]         # A or B
        date = parts[1]           # 2022-09-01
        replicate = parts[-1]     # 1 or 2

        return prefix, date, replicate

    def is_valid_pair(self, row):
        pre = str(row["preselection_sample"]).strip()
        post = str(row["postselection_sample"]).strip()

        # ✔ Only keep correct experiment types
        if "VSVG_control" not in pre:
            return False
        if "no-antibody_control" not in post:
            return False

        # ✔ Extract metadata
        try:
            p_prefix, p_date, p_rep = self.parse_sample(pre)
            q_prefix, q_date, q_rep = self.parse_sample(post)
        except Exception:
            return False

        #  Enforce strict matching
        return (
            p_prefix == q_prefix and
            p_date == q_date and
            p_rep == q_rep
        )

    def clean(self):
        print("\nLoading functional selections...")
        df = pd.read_csv(self.input_csv)

        print(f"Total original pairs: {len(df)}")

        # Apply filtering
        clean_df = df[df.apply(self.is_valid_pair, axis=1)].copy()

        # Remove duplicates (extra safety)
        clean_df = clean_df.drop_duplicates(
            subset=["preselection_sample", "postselection_sample"]
        )

        print(f"Total valid pairs after cleaning: {len(clean_df)}")

        # Prevent overwrite
        if os.path.exists(self.output_csv):
            raise FileExistsError(
                f"Output file already exists: {self.output_csv}\n"
                "Delete it or change the name."
            )

        # Save
        clean_df.to_csv(self.output_csv, index=False)

        print("\n Clean file saved successfully:")
        print(self.output_csv)


#  RUN SCRIPT
#if __name__ == "__main__":

    #pairer = CleanPairer(
        #input_csv="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/func_scores/functional_selections.csv",
        #output_csv="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/func_scores/functional_selections_clean.csv"
    #)

    #pairer.clean()