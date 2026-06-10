import pandas as pd
import os


class VariantMerger:
    def __init__(self, variant_count_dir, selection_file, output_folder="merged_output"):

        self.variant_count_dir = variant_count_dir
        self.selection_file = selection_file

        # Portable output path (based on selection file location)
        base_dir = os.path.dirname(self.selection_file)
        self.output_dir = os.path.join(base_dir, output_folder)

        os.makedirs(self.output_dir, exist_ok=True)

        self.selections = pd.read_csv(self.selection_file)

        print("Variant count dir:", self.variant_count_dir)
        print("Selection file:", self.selection_file)
        print("Output dir:", self.output_dir)

    def load_and_filter(self, filepath):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"File not found: {filepath}")

        df = pd.read_csv(filepath)


        #check that  mising columns in a csv file
        required_columns = ['barcode', 'count', 'variant_call_support', 'aa_substitutions', 'codon_substitutions']

        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            raise ValueError(f"Missing columns in {filepath}: {missing_columns}")




        # Remove neut_standard rows
        df = df[
            (df["aa_substitutions"] != "neut_standard") &
            (df["codon_substitutions"] != "neut_standard")
        ]


        return df
    


    def find_file(self, sample_name, selection_name):
        prefix = selection_name.split("_")[0]
        matches = [
            f for f in os.listdir(self.variant_count_dir)
            if f.endswith(".csv") 
            and f.startswith(prefix + "_") 
            and sample_name in f
        ]

        if len(matches) == 0:
            raise FileNotFoundError(f"No file found for sample: {sample_name}")

        if len(matches) > 1:
            raise ValueError(
            f"Multiple matches found for {sample_name}: {matches}\n"
            " You need to refine matching (likely A_ vs B_ issue)"
        )

        return os.path.join(self.variant_count_dir, matches[0])




    def merge_variant_counts(self, pre_df, post_df):
        # Rename columns to avoid collisions
        pre_df = pre_df.rename(columns={
            "count": "pre_count",
            "variant_call_support": "pre_variant_call_support"
        })

        post_df = post_df.rename(columns={
            "count": "post_count",
            "variant_call_support": "post_variant_call_support"
        })


        # drop duplicates
        pre_df = pre_df.drop_duplicates(subset=["barcode", "aa_substitutions", "codon_substitutions"])
        post_df = post_df.drop_duplicates(subset=["barcode", "aa_substitutions", "codon_substitutions"])



        # Merge on unique identifiers
        merged = pd.merge(
            pre_df,
            post_df,
            on=["barcode", "aa_substitutions", "codon_substitutions"],
            how="outer"  # change to "outer" if needed
        )



        # Fill missing values with 0
        merged["pre_count"] = merged["pre_count"].fillna(0)
        merged["post_count"] = merged["post_count"].fillna(0)


        #fill missing variant call support with 0

        merged["pre_variant_call_support"] = merged["pre_variant_call_support"].fillna(0)
        merged["post_variant_call_support"] = merged["post_variant_call_support"].fillna(0)



        #convert variant call support to integers
        merged["pre_variant_call_support"] = merged["pre_variant_call_support"].astype(int)
        merged["post_variant_call_support"] = merged["post_variant_call_support"].astype(int)   



        # Convert counts back to integers
        merged["pre_count"] = merged["pre_count"].astype(int)
        merged["post_count"] = merged["post_count"].astype(int)

        return merged

    def process_single_selection(self, row):

        selection_name = str(row["selection_name"]).strip()
        pre_file = str(row["preselection_sample"]).strip()
        post_file = str(row["postselection_sample"]).strip()




        pre_path = self.find_file(pre_file, selection_name)
        post_path = self.find_file(post_file, selection_name)   



        print(f"\nProcessing: {selection_name}")
        print("Pre path:", pre_path)
        print("Post path:", post_path)

        # Load and filter
        pre_df = self.load_and_filter(pre_path)
        post_df = self.load_and_filter(post_path)

        # Merge
        merged = self.merge_variant_counts(pre_df, post_df)

        # Select required columns
        final_df = merged[[
            "barcode",
            "pre_count",
            "post_count",
            "codon_substitutions",
            "aa_substitutions",
            "pre_variant_call_support",
            "post_variant_call_support"
        ]].copy()  # avoid pandas warning

        # Add selection name
        final_df["selection_name"] = selection_name

        # Reorder columns
        final_df = final_df[[
            "selection_name",
            "barcode",
            "pre_count",
            "post_count",
            "codon_substitutions",
            "aa_substitutions",
            "pre_variant_call_support",
            "post_variant_call_support"
        ]]

        return selection_name, final_df

    def run_all(self):
        results = {}

        for _, row in self.selections.iterrows():
            selection_name, df = self.process_single_selection(row)

            # Save each selection as its own CSV
            safe_name = selection_name.replace(" ", "_")
            output_path = os.path.join(self.output_dir, f"{safe_name}_merged.csv")
            df.to_csv(output_path, index=False)

            results[selection_name] = df

        print("\n All selections processed successfully.")
        print(f" Files saved in: {self.output_dir}")

        return results


#  OUTSIDE CLASS (NO INDENTATION)
#if __name__ == "__main__":

    #merger = VariantMerger(
        #variant_count_dir="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/variant_counts",
        #selection_file="/home/lechiffre/HIV_Envelope_BF520_DMS_CD4bs_sera/results/func_scores/functional_selections_clean.csv",
        #output_folder="merged_output_clean"
    #)

    #merger.run_all()