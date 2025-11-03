#!/usr/bin/env python
import argparse
import re

def get_umi(header: str) -> str:
    """
    Takes a SAM read header and returns the corresponding UMI.
    """
    if not header:
        return ""
    
    parts = header.split(':')
    return parts[-1]

def get_pos_updated(position_orig: int, forward_stranded: bool, cigar: str) -> int:
    """
    Takes the original position, strand, and CIGAR string and outputs the
    adjusted 5' position index, accounting for soft clippings. 

    1.  Input: 10, true, 3S5M2S -> Output: 7
        (Forward strand: pos - leading_soft_clip = 10 - 3 = 7)
    2.  Input: 10, false, 3S5M2S -> Output: 17
        (Reverse strand: pos + ref_span + trailing_soft_clip = 10 + 5 + 2 = 17)
    """
    
    # Use regex to parse CIGAR string 
    cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    
    if not cigar_ops:
        return position_orig

    if forward_stranded:
        #Check the first operation
        first_op_val_str, first_op_code = cigar_ops[0]
        if first_op_code == 'S':
            #Subtract leading soft clips
            return position_orig - int(first_op_val_str)
        else:
            #No leading soft clips
            return position_orig
    else:
        #Reverse stranded
        pos_updated = position_orig
        
        #Operations that consume the reference sequence
        ref_consuming_ops = {'M', 'D', 'N', 'X', '='}
        
        #Add the length of all reference-consuming operations
        for val_str, code in cigar_ops:
            if code in ref_consuming_ops:
                pos_updated += int(val_str)
        
        #Check the last operation for trailing soft clips
        last_op_val_str, last_op_code = cigar_ops[-1]
        if last_op_code == 'S':
            #Add trailing soft clips
            pos_updated += int(last_op_val_str)
            
        return pos_updated

def get_args():
    parser = argparse.ArgumentParser(
        description=(
            "Usage: israel_deduper.py [-h] [-f STRING] [-o STRING] [-u STRING]\n"
            "Given a SAM file of uniquely mapped reads, remove all PCR duplicates (retain only the first detected single copy of each read):\n"
            "NOTE: This script assumes SAM file contains unique paired end reads, UMI length of 8, and uncompressed file inputs.\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter
        )
    parser.add_argument(
        "-f", "--file",
        required=True,
        help="designates absolute file path to sorted sam file"
    )
    parser.add_argument(
        "-o", "--outfile",
        required=True,
        help="designates absolute file path to deduplicated sam file"
    )
    parser.add_argument(
        "-u", "--umi",
        required=True,
        help="designates file containing the list of UMIs"
    )
    return parser.parse_args()

def main():
    args = get_args() 

    umi_set = set()
    with open(args.umi, 'r') as f_umi: # Input UMI file and make a set of UMIs.
        for umi in f_umi:
            umi_set.add(umi.strip())
    
    current_rname = None
    seen_keys = set()

    # Statistics counters
    header_count = 0
    total_alignments = 0
    wrong_umi_count = 0
    duplicate_count = 0
    unique_reads_count = 0
    chrom_counts = {}

    with open(args.file, 'r') as sam_in, open(args.outfile, 'w') as sam_out:
        for line in sam_in:
            if line.startswith('@'): # If line is a header line, write to output automatically
                header_count += 1
                sam_out.write(line)
                continue
            
            total_alignments += 1

            fields = line.strip().split('\t')
   
            if len(fields) < 11: # Basic check for a valid SAM alignment line
                continue
            
            header = fields[0]
            flag = int(fields[1])
            rname = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            umi = get_umi(header)
            forward_stranded = (flag & 16) != 16

            if umi not in umi_set: # Unknown umi, skip
                wrong_umi_count += 1
                continue
            
            if current_rname is None: # First chromosome detected
                current_rname = rname
            elif rname != current_rname: # New chromosome detected, clear memory for previous chromosome
                seen_keys.clear()
                current_rname = rname
                
            pos_updated = get_pos_updated(pos, forward_stranded, cigar) # Get updated 5' position
            k = (rname, pos_updated, forward_stranded, umi)

            if k not in seen_keys:
                unique_reads_count += 1
                sam_out.write(line)
                seen_keys.add(k)
                chrom_counts[rname] = chrom_counts.get(rname, 0) + 1
            else:
                duplicate_count += 1

        print("\n--- Summary ---")
        print(f"Header lines: \t{header_count}")
        print(f"Total alignments processed: \t{total_alignments}")
        print(f"Reads with wrong UMIs: \t{wrong_umi_count}")
        print(f"Duplicates removed: \t{duplicate_count}")
        print(f"Unique reads written: \t{unique_reads_count}")

        print("\n--- Unique Reads Per Chromosome ---")
        for chrom, count in chrom_counts.items():
            print(f"{chrom}\t{count}")

if __name__ == "__main__":
    main()