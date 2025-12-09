import pandas as pd

def read_vcf(vcf_file):
    vcf_data = []

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) >= 5:
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]

                var_type = "SNP" if len(ref) == len(alt) else "INDEL"

                vcf_data.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "Type": var_type,
                    "Source": "bcftools"
                })

    return pd.DataFrame(vcf_data)

def read_csv(csv_file):
    df = pd.read_csv(csv_file)

    df_renamed = df.rename(columns={
        "POS": "POS",
        "REF": "REF",
        "ALT": "ALT"
    })

    df_renamed["POS"] = df_renamed["POS"].astype(int)
    df_renamed["Type"] = df_renamed.apply(
        lambda row: "SNP" if len(str(row["REF"])) == len(str(row["ALT"])) else "INDEL",
        axis=1
    )
    df_renamed["Source"] = "CSV"

    return df_renamed[["POS", "REF", "ALT", "Type", "Source"]]

def merge_and_compare(vcf_file, csv_file, output_file):

    vcf_df = read_vcf(vcf_file)
    csv_df = read_csv(csv_file)

    vcf_df["variant_key"] = vcf_df["POS"].astype(str) + "_" + vcf_df["REF"] + "_" + vcf_df["ALT"]
    csv_df["variant_key"] = csv_df["POS"].astype(str) + "_" + csv_df["REF"] + "_" + csv_df["ALT"]

    vcf_keys = set(vcf_df["variant_key"])
    csv_keys = set(csv_df["variant_key"])

    matched_keys = vcf_keys & csv_keys
    vcf_only_keys = vcf_keys - csv_keys
    csv_only_keys = csv_keys - vcf_keys

    result_rows = []

    for key in matched_keys:
        row = vcf_df[vcf_df["variant_key"] == key].iloc[0]
        result_rows.append({
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "Type": row["Type"],
            "Match_Status": "MATCH",
            "Source": "Both"
        })

    for key in vcf_only_keys:
        row = vcf_df[vcf_df["variant_key"] == key].iloc[0]
        result_rows.append({
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "Type": row["Type"],
            "Match_Status": "MISMATCH",
            "Source": "VCF_only"
        })

    for key in csv_only_keys:
        row = csv_df[csv_df["variant_key"] == key].iloc[0]
        result_rows.append({
            "POS": row["POS"],
            "REF": row["REF"],
            "ALT": row["ALT"],
            "Type": row["Type"],
            "Match_Status": "MISMATCH",
            "Source": "CSV_only"
        })

    result_df = pd.DataFrame(result_rows)
    result_df = result_df.sort_values("POS").reset_index(drop=True)
    result_df.to_csv(output_file, index=False)

    return result_df

if __name__ == "__main__":
    vcf_file = "variants.vcf"
    csv_file = "simulated_mutated_genome.csv"
    output_file = "merged_result.csv"

    result = merge_and_compare(vcf_file, csv_file, output_file)
