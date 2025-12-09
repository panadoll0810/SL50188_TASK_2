import random

def genome_to_paired_reads(genome_file, output_prefix, read_length=100, num_reads=None,
                           coverage=30, quality_score="I", insert_size=300):
    """prepare paired-end reads based on simulated mutated genome"""
    genome_seq = read_genome(genome_file)

    if not genome_seq:
        print("Error: Can NOT read fasta file")
        return

    genome_length = len(genome_seq)

    if num_reads is None:
        num_reads = int((genome_length * coverage) / (read_length * 2))

    reads_r1 = []
    reads_r2 = []
    read_count = 0

    while read_count < num_reads:
        max_start = max(0, genome_length - insert_size - read_length)
        if max_start <= 0:
            break

        start_pos = random.randint(0, max_start)

        read1_seq = genome_seq[start_pos:start_pos + read_length]
        if len(read1_seq) != read_length:
            continue

        read2_start = start_pos + insert_size
        read2_seq = genome_seq[read2_start:read2_start + read_length]
        if len(read2_seq) != read_length:
            continue

        read2_seq = reverse_complement(read2_seq)

        read_id = f"read_{read_count + 1}"

        reads_r1.append({
            "id": f"{read_id}/1",
            "sequence": read1_seq,
            "quality": quality_score * len(read1_seq)
        })

        reads_r2.append({
            "id": f"{read_id}/2",
            "sequence": read2_seq,
            "quality": quality_score * len(read2_seq)
        })

        read_count += 1

    r1_file = f"{output_prefix}_R1.fastq"
    r2_file = f"{output_prefix}_R2.fastq"

    write_fastq(r1_file, reads_r1)
    write_fastq(r2_file, reads_r2)

    print(f"Saved {r1_file}")
    print(f"Saved {r2_file}")


def read_genome(genome_file):
    """read genome file
        can be used on both txt and fasta files"""
    sequence = []

    try:
        with open(genome_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(">"):
                    continue
                sequence.append(line)

        if not sequence:
            print("Error：Can NOT Find Sequence")
            return None

        return "".join(sequence).upper()

    except FileNotFoundError:
        print(f"Error：Can NOT Find File {genome_file}")
        return None

    except Exception as e:
        print(f"Error：Reading Files - {e}")
        return None

def reverse_complement(seq):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq))

def write_fastq(output_file, reads):
    with open(output_file, "w") as f:
        for read in reads:
            if len(read["sequence"]) != len(read["quality"]):
                print(f"Warning: {read['id']} Sequence Length and Quality Score Length Do NOT Match")
                continue

            f.write(f"@{read['id']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            f.write(f"{read['quality']}\n")

if __name__ == "__main__":
    genome_to_paired_reads(
        "/Users/milkcaramelcheng/Desktop/task2/E_coli/simulated_mutate_genome.txt",
        "/Users/milkcaramelcheng/Desktop/task2/E_coli/simulated_read",
        read_length=100,
        coverage=30,
        insert_size=300
    )