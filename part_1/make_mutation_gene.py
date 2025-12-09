import random
import re
import csv

class SafeStringEditor:
    def __init__(self, original_string):
        self.string = original_string
        self.original_string = original_string
        
        # insert markers to wrap insertion
        self.marker_prefix = "<<<IN_"
        self.marker_suffix = "_IN>>>"
        self.insert_counter = 0
        
        # protect SNPs from deletion
        self.protected_chars = {"@", "#", ".", "%"}
        self.snp_map = {}
        self.snp_records = []
        self.operations_log = []
        self.position_mapping = list(range(len(original_string)))

    def pre_snps_sequence(self, count):
        
        """randomly label the positions of SNP by replace them with special characters
            protect them from deletion and avoid any changes of the insertion"""

        if not self.string:
            return self

        seq = list(self.string)
        count = min(count, len(seq))
        positions = random.sample(range(len(seq)), count)

        replace_char = {
            "A": "@",
            "T": "#",
            "G": ".",
            "C": "%"
        }

        for pos in positions:
            char = seq[pos].upper()
            if char in replace_char:
                self.snp_map[pos] = {
                    "original": char,
                    "symbol": replace_char[char],
                    "candidates": [c for c in "ATGC" if c != char]
                }
                seq[pos] = replace_char[char]

        self.string = "".join(seq)
        return self

    def get_marker_positions(self):
        
        """get positions of all insertion markers"""
        
        pattern = re.escape(self.marker_prefix) + r"\d+" + re.escape(self.marker_suffix) + r".*?" + re.escape(
            self.marker_prefix) + r"\d+_END" + re.escape(self.marker_suffix)
        insertions = list(re.finditer(pattern, self.string))

        marker_positions = set()
        for match in insertions:
            for pos in range(match.start(), match.end()):
                marker_positions.add(pos)
        return marker_positions

    def string_pos_to_mapping_index(self, string_pos):
       
        """convert current string position to index in position mapping table
            skip marker positions"""
        
        marker_positions = self.get_marker_positions()
        mapping_index = 0
        for i in range(string_pos):
            if i not in marker_positions:
                mapping_index += 1
        return mapping_index

    def restore_snps(self):
        
        """Restore labelled SNPs
            replace the special characters with random selected bases"""
        
        replace_back_char = {
            "@": ["T", "G", "C"],
            "#": ["A", "G", "C"],
            ".": ["A", "T", "C"],
            "%": ["A", "T", "G"],
        }

        marker_positions = self.get_marker_positions()
        seq = list(self.string)
        self.snp_records = []

        mapping_index = 0
        for current_pos, char in enumerate(seq):
            if current_pos in marker_positions:
                continue

            if char in replace_back_char:
                chosen = random.choice(replace_back_char[char])
                seq[current_pos] = chosen

                if mapping_index < len(self.position_mapping):
                    original_pos = self.position_mapping[mapping_index]

                    if original_pos in self.snp_map:
                        self.snp_records.append({
                            "position": original_pos,
                            "ref": self.snp_map[original_pos]["original"],
                            "alt": chosen
                        })

            mapping_index += 1

        self.string = "".join(seq)
        return self

    def get_safe_insertion_ranges(self):
        
        """get list of positions where insertions are safe"""
        
        pattern = re.escape(self.marker_prefix) + r"\d+" + re.escape(self.marker_suffix) + r".*?" + re.escape(
            self.marker_prefix) + r"\d+_END" + re.escape(self.marker_suffix)
        insertions = list(re.finditer(pattern, self.string))

        safe_positions = []
        last_end = 0

        for match in insertions:
            safe_positions.extend(range(last_end, match.start() + 1))
            last_end = match.end()

        safe_positions.extend(range(last_end, len(self.string) + 1))
        return safe_positions

    def get_actual_position(self, string_pos):
        
        """get the actual position in original sequence for current string position"""
        
        mapping_index = self.string_pos_to_mapping_index(string_pos)
        if mapping_index < len(self.position_mapping):
            return self.position_mapping[mapping_index]
        return mapping_index

    def insert_random(self, sequence):
       
        """insert a sequence at a random position"""
        
        safe_positions = self.get_safe_insertion_ranges()

        pos = random.choice(safe_positions)
        marked_seq = f"{self.marker_prefix}{self.insert_counter}{self.marker_suffix}{sequence}{self.marker_prefix}{self.insert_counter}_END{self.marker_suffix}"
        self.string = self.string[:pos] + marked_seq + self.string[pos:]
        actual_pos = self.get_actual_position(pos)

        self.operations_log.append({
            "type": "insert",
            "position": actual_pos,
            "sequence": sequence,
            "length": len(sequence)
        })

        self.insert_counter += 1
        return actual_pos

    def get_safe_deletion_ranges(self):
        
        """get list of ranges where deletions are safe"""
        
        marker_positions = self.get_marker_positions()
        unsafe_positions = set(marker_positions)

        for i, char in enumerate(self.string):
            if char in self.protected_chars:
                unsafe_positions.add(i)

        safe_ranges = []
        current_start = None

        for i in range(len(self.string)):
            if i not in unsafe_positions:
                if current_start is None:
                    current_start = i
            else:
                if current_start is not None:
                    safe_ranges.append((current_start, i))
                    current_start = None

        if current_start is not None:
            safe_ranges.append((current_start, len(self.string)))

        return safe_ranges

    def delete_random(self, length):
        
        """delete a random sequence"""
        
        safe_ranges = self.get_safe_deletion_ranges()

        total_safe_length = sum(end - start for start, end in safe_ranges)

        if total_safe_length < length:
            print(f"Warning: Only {total_safe_length} chars could be deleted, requested {length}")
            length = total_safe_length

        if total_safe_length == 0:
            print("Warning: No chars could be deleted")
            return None

        del_start = random.randint(0, total_safe_length - length)

        accumulated = 0
        for start, end in safe_ranges:
            range_len = end - start
            if accumulated + range_len > del_start:
                offset = del_start - accumulated
                actual_start = start + offset
                actual_end = min(actual_start + length, end)

                deleted_chars = self.string[actual_start:actual_end]
                actual_pos = self.get_actual_position(actual_start)

                start_mapping_index = self.string_pos_to_mapping_index(actual_start)
                end_mapping_index = self.string_pos_to_mapping_index(actual_end)

                del self.position_mapping[start_mapping_index:end_mapping_index]

                self.string = self.string[:actual_start] + self.string[actual_end:]

                self.operations_log.append({
                    "type": "delete",
                    "position": actual_pos,
                    "deleted_chars": deleted_chars,
                    "length": len(deleted_chars)
                })
                return actual_pos
            accumulated += range_len

        return None

    def perform_indels(self, num_indels):
        
        """perform specified number of insertion and deletion"""
        
        operations = []
        for i in range(num_indels):
            if random.choice(["insert", "delete"]) == "insert":
                insert_length = random.randint(1, 10)
                sequence = "".join(random.choices("ATGC", k=insert_length))
                operations.append(("insert", sequence))
            else:
                delete_length = random.randint(1, 10)
                operations.append(("delete", delete_length))

        for op_type, param in reversed(operations):
            if op_type == "insert":
                self.insert_random(param)
            else:
                self.delete_random(param)

        return

    def get_final_string(self):
        
        """get final string
            remove all markers"""
        
        pattern = re.escape(self.marker_prefix) + r"\d+" + re.escape(self.marker_suffix)
        result = re.sub(pattern, "", self.string)
        pattern_end = re.escape(self.marker_prefix) + r"\d+_END" + re.escape(self.marker_suffix)
        result = re.sub(pattern_end, "", result)
        return result

    def generate_report(self):
        
        """generate a report of the mutations
            record positions of mutations
            record types of mutations"""
        
        report = []
        report.append(f"{'Operation':<12} {'POS':<10} {'REF':<50} {'ALT':<50}")

        all_operations = []

        for snp in self.snp_records:
            all_operations.append({
                "type": "SNP",
                "position": snp["position"] + 1,  # 1-based position
                "ref": snp["ref"],
                "alt": snp["alt"]
            })

        for op in self.operations_log:
            op_type = op["type"]
            position = op["position"]

            if op_type == "insert":
                if position > 0:
                    prev_char = self.original_string[position - 1]
                    ref = prev_char
                    alt = prev_char + op["sequence"]
                else:
                    ref = ""
                    alt = op["sequence"]
                op_display = "Insertion"
            else:
                if position > 0:
                    prev_char = self.original_string[position - 1]
                    ref = prev_char + op["deleted_chars"]
                    alt = prev_char
                else:
                    ref = op["deleted_chars"]
                    alt = ""
                op_display = "Deletion"

            all_operations.append({
                "type": op_display,
                "position": position,
                "ref": ref,
                "alt": alt
            })

        all_operations.sort(key=lambda x: x["position"])

        for op in all_operations:
            ref_display = op["ref"] if len(op["ref"]) <= 50 else op["ref"][:47] + "..."
            alt_display = op["alt"] if len(op["alt"]) <= 50 else op["alt"][:47] + "..."
            report.append(f"{op['type']:<12} {op['position']:<10} {ref_display:<50} {alt_display:<50}")

        with open("simulated_mutated_genome.txt", "w") as f:
            f.write(self.get_final_string())

        return "\n".join(report)

def read_fasta_sequence(fasta_file):
    try:
        with open(fasta_file, "r") as f:
            lines = f.readlines()
            if len(lines) < 2:
                print(f"Error: File {fasta_file} has insufficient content")
                return None
            sequence = "".join([line.strip() for line in lines[1:]])
            return sequence
    except FileNotFoundError:
        print(f"Error: File {fasta_file} not found")
        return None

if __name__ == "__main__":
    original = read_fasta_sequence("reference_genome.fasta")
    editor = SafeStringEditor(original)
    editor.pre_snps_sequence(300)
    editor.perform_indels(20)
    editor.restore_snps()

    report = editor.generate_report()
    print(report)

    """store the mutations in a csv file"""
    with open("simulated_mutated_genome.csv", "w", newline="") as f:
        lines = [line.split() for line in report.split("\n") if line]
        header = lines[0] if lines else []
        data = lines[1:] if len(lines) > 1 else []

        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)

