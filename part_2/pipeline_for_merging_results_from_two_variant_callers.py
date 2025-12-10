import subprocess
import sys
import pandas as pd
import pysam


def run_command(cmd_list, description):
    print(f"Running: {description}")

    try:
        result = subprocess.run(cmd_list, check=True)
        print(f"{description} Completed\n")
        return True

    except subprocess.CalledProcessError as e:
        print(f"Error in {description}: {e}", file=sys.stderr)
        return False

    except FileNotFoundError as e:
        print(f"Command NOT Found: {e}", file=sys.stderr)
        return False


def run_pipeline(proc1_cmd, proc2_cmd, proc3_cmd, output_file, description):
    print(f"Running: {description}")

    try:
        proc1 = subprocess.Popen(proc1_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        proc2 = subprocess.Popen(
            proc2_cmd,
            stdin=proc1.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        proc1.stdout.close()

        with open(output_file, "wb") as f_out:
            proc3 = subprocess.Popen(
                proc3_cmd,
                stdin=proc2.stdout,
                stdout=f_out,
                stderr=subprocess.PIPE
            )
            proc2.stdout.close()

            proc3.wait()
            proc2.wait()
            proc1.wait()

        if proc1.returncode != 0:
            print(f"Process 1 Error: {proc1.stderr.read().decode()}", file=sys.stderr)
            return False

        if proc2.returncode != 0:
            print(f"Process 2 Error: {proc2.stderr.read().decode()}", file=sys.stderr)
            return False

        if proc3.returncode != 0:
            print(f"Process 3 Error: {proc3.stderr.read().decode()}", file=sys.stderr)
            return False

        print(f"{description} Completed successfully")
        return True

    except Exception as e:
        print(f"Pipeline Error: {e}", file=sys.stderr)
        return False


def main():
    # File paths
    reference = "Ecoli_complete_genome.fasta"
    reference_bcf = "Ecoli_complete_genome.fasta"
    read_r1 = "SRR25083113_1.fastq.gz"
    read_r2 = "SRR25083113_2.fastq.gz"
    bam_file = "Ecoli_mapped_reads.bam"
    vcf_bcf = "Ecoli_variants_bcf.vcf"
    vcf_bcf_gz = "Ecoli_variants_bcf.vcf.gz"
    snippy_dir = "Ecoli_snippy_results"
    snippy_vcf_gz = "Ecoli_snippy_results/snps.vcf.gz"
    merged_vcf = "Ecoli_merged.vcf.gz"

    # Step 1: Map reads (minimap2 | samtools view | samtools sort)
    map_pipeline = run_pipeline(
        ["minimap2", "-a", "-x", "sr", reference, read_r1, read_r2],
        ["samtools", "view", "-h", "-F", "0x900", "-"],
        ["samtools", "sort", "-O", "bam"],
        bam_file,
        "Map reads with minimap2, and convert/sort the output into a BAM file using samtools."
    )
    if not map_pipeline:
        sys.exit(1)

    # Step 2: Call variants with bcftools
    print("Running: Call variants with bcftools")

    try:
        mpileup_proc = subprocess.Popen(
            ["bcftools", "mpileup", "-Ou", "-f", reference_bcf, bam_file],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        with open(vcf_bcf, 'w') as vcf_out:
            call_proc = subprocess.Popen(
                ["bcftools", "call", "-vc", "-Ov"],
                stdin=mpileup_proc.stdout,
                stdout=vcf_out,
                stderr=subprocess.PIPE
            )
            mpileup_proc.stdout.close()

            mpileup_stderr = mpileup_proc.communicate()[1].decode()
            call_stdout, call_stderr = call_proc.communicate()

        if mpileup_proc.returncode != 0:
            print(f"bcftools mpileup Error:\n{mpileup_stderr}", file=sys.stderr)
            sys.exit(1)

        if call_proc.returncode != 0:
            print(f"bcftools call Error:\n{call_stderr.decode()}", file=sys.stderr)
            sys.exit(1)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Run snippy
    if not run_command(
            ["snippy", "--outdir", snippy_dir, "--ref", reference,
             "--R1", read_r1, "--R2", read_r2, "--mincov", "10", "--minfrac", "0.9"],
            "Run snippy variant calling"
    ):
        sys.exit(1)

    # Step 4: Compress VCF with bgzip
    if not run_command(["bgzip", vcf_bcf], "Compress VCF with bgzip"):
        sys.exit(1)

    # Step 5: Create index for bcftools VCF
    if not run_command(
            ["tabix", "-p", "vcf", vcf_bcf_gz],
            "Create tabix index for bcftools VCF"
    ):
        sys.exit(1)

    # Step 6: Create index for snippy VCF
    if not run_command(
            ["tabix", "-p", "vcf", snippy_vcf_gz],
            "Create tabix index for snippy VCF"
    ):
        sys.exit(1)

    # Step 7: Merge VCF files
    if not run_command(
            ["bcftools", "merge", vcf_bcf_gz, snippy_vcf_gz, "-Oz", "-o", merged_vcf],
            "Merge VCF files with bcftools"
    ):
        sys.exit(1)

    # Step 8: Extract genotypes from merged VCF
    try:
        vcf = pysam.VariantFile(merged_vcf)
        records = []
        samples = list(vcf.header.samples)

        def get_allele(gt, rec):

            if not gt or gt[0] is None:
                return "missing"

            allele_indices = []
            if isinstance(gt, (tuple, list)):
                for allele in gt:
                    if allele is not None:
                        allele_indices.append(allele)
            else:
                allele_indices.append(gt)

            if not allele_indices:
                return "missing"

            non_ref = []
            for idx in allele_indices:
                if idx > 0:
                    non_ref.append(idx)

            if non_ref:
                first_mutation_idx = non_ref[0]
                return rec.alts[first_mutation_idx - 1]
            else:
                return rec.ref

        for rec in vcf:
            row = [rec.chrom, rec.pos, rec.ref, ",".join(rec.alts)]

            for sample in samples:
                gt = rec.samples[sample].get("GT")
                allele = get_allele(gt, rec)
                row.append(allele)
            records.append(row)
        vcf.close()

        columns = ["CHROM", "POS", "REF", "ALT"]
        columns.extend(samples)
        df = pd.DataFrame(records, columns=columns)
        df.to_csv("merged_Ecoli.csv", index=False)

    except Exception as e:
        print(f"Error extracting genotypes: {e}", file=sys.stderr)
        sys.exit(1)

    print("Complete pipeline finished")


if __name__ == "__main__":
    main()