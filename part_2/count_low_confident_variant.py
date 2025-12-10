import gzip

input_vcf = "Ecoli_merged.vcf.gz"
output_csv = "not_confident_calls.csv"

with gzip.open(input_vcf, "rt") as f_in, open(output_csv, "w") as f_out:
    header = ["CHROM", "POS", "REF", "ALT", "QUAL", "GT_bcf", "GT_snippy", "Reason"]
    f_out.write(",".join(header) + "\n")

    count = 0

    for line in f_in:
        if line.startswith("#"):
            continue

        cols = line.strip().split("\t")

        chrom = cols[0]
        pos = cols[1]
        ref = cols[3]
        alt = cols[4]
        qual_str = cols[5]

        try:
            qual = float(qual_str)
        except ValueError:
            qual = 0.0

        format_parts = cols[8].split(":")
        try:
            gt_index = format_parts.index("GT")
        except ValueError:
            continue

        sample1_data = cols[9].split(":")
        sample2_data = cols[10].split(":")

        gt1 = sample1_data[gt_index]
        gt2 = sample2_data[gt_index]

        reasons = []

        if qual < 20:
            reasons.append("Low_QUAL")

        if gt1 == "0/1" or gt2 == "0/1":
            reasons.append("Fake_Heterozygous")

        if gt1 != gt2:
            reasons.append("Tool_Mismatch")

        if len(reasons) > 0:
            count += 1
            row = [
                chrom,
                pos,
                ref,
                alt,
                str(qual),
                f"'{gt1}",
                f"'{gt2}",
                ";".join(reasons)
            ]
            f_out.write(",".join(row) + "\n")

print(f"found {count} non_confident variant")
print(f"saved to {output_csv}")