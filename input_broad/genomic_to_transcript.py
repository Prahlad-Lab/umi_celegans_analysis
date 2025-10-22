import argparse
import sys
from collections import defaultdict

import pysam
from intervaltree import Interval, IntervalTree

def parse_gtf(gtf_file_path):
    """
    Parses a GTF file to build transcript structures and an interval tree for fast lookups.

    Returns:
        tuple: Contains (transcript_map, interval_trees, all_transcript_ids)
    """
    print("Parsing GTF file to build transcript map and interval trees...")
    transcripts = defaultdict(lambda: {'exons': [], 'strand': None})
    interval_trees = defaultdict(IntervalTree)
    all_transcript_ids = set()

    with open(gtf_file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 8 or fields[2] != 'exon':
                continue

            chrom, _, _, start, end, _, strand, _, attributes = fields
            start, end = int(start), int(end)

            match = re.search(r'transcript_id "([^"]+)"', attributes)
            if not match:
                continue
            transcript_id = match.group(1)
            
            all_transcript_ids.add(transcript_id)
            transcripts[transcript_id]['strand'] = strand
            transcripts[transcript_id]['exons'].append((start, end))
            # Store transcript_id in the interval tree for easy retrieval
            interval_trees[chrom].add(Interval(start, end + 1, transcript_id))

    # Sort exons by genomic position for each transcript
    for tid in transcripts:
        transcripts[tid]['exons'].sort(key=lambda x: x[0])

    print(f"Parsed {len(transcripts)} transcripts.")
    return dict(transcripts), interval_trees, all_transcript_ids

def reverse_complement(seq):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement.get(base, base) for base in reversed(seq))

def convert_vcf_coordinates(vcf_path, gtf_path, output_path, prefix):
    """
    Removes a prefix from chromosome names and converts genomic coordinates to transcriptomic.
    """
    transcripts, interval_trees, all_transcript_ids = parse_gtf(gtf_path)

    try:
        vcf_in = pysam.VariantFile(vcf_path, "r")
    except (IOError, ValueError) as e:
        print(f"Error opening input VCF file: {e}", file=sys.stderr)
        sys.exit(1)

    # Build a new header with transcript IDs as contigs
    new_header = pysam.VariantHeader()
    for record in vcf_in.header.records:
        if record.key != 'contig':
            new_header.add_record(record)

    for tid in sorted(list(all_transcript_ids)):
        new_header.add_line(f"##contig=<ID={tid}>")
    
    new_header.info.add("OriginalGenomicChrom", "1", "String", "Original genomic chromosome of the variant.")
    new_header.info.add("OriginalGenomicPos", "1", "Integer", "Original genomic position of the variant.")

    for sample in vcf_in.header.samples:
        new_header.add_sample(sample)

    vcf_out = pysam.VariantFile(output_path, "w", header=new_header)
    
    print(f"Starting conversion. Removing prefix '{prefix}' from chromosome names...")
    converted_count = 0
    skipped_count = 0

    for record in vcf_in:
        # 1. Modify chromosome name
        original_chrom = record.chrom
        genomic_chrom = original_chrom.replace(prefix, '')
        genomic_pos = record.pos

        # 2. Find overlapping transcripts using the interval tree
        if genomic_chrom not in interval_trees:
            skipped_count += 1
            continue
            
        overlapping_exons = interval_trees[genomic_chrom].at(genomic_pos)
        if not overlapping_exons:
            # Variant is in an intron or intergenic region
            skipped_count += 1
            continue

        # 3. For each mapping, calculate the new coordinates
        for exon_interval in overlapping_exons:
            transcript_id = exon_interval.data
            tx_data = transcripts[transcript_id]
            
            transcriptomic_pos = 0
            strand = tx_data['strand']
            
            # Use a copy of exons sorted appropriately for transcriptomic calculation
            exons_for_calc = tx_data['exons']
            if strand == '-':
                # For negative strand, reverse the exon order for calculation
                exons_for_calc = sorted(tx_data['exons'], key=lambda x: x[0], reverse=True)

            found = False
            for exon_start, exon_end in exons_for_calc:
                exon_len = exon_end - exon_start + 1
                if exon_start <= genomic_pos <= exon_end:
                    # This is the correct exon
                    if strand == '+':
                        offset = genomic_pos - exon_start
                        transcriptomic_pos += offset + 1
                    else: # strand == '-'
                        offset = exon_end - genomic_pos
                        transcriptomic_pos += offset + 1
                    found = True
                    break
                else:
                    transcriptomic_pos += exon_len
            
            if not found: continue # Should not happen if logic is correct

            # 4. Write the new record
            new_record = vcf_out.new_record(contig=transcript_id, start=transcriptomic_pos - 1)
            new_record.id = record.id
            new_record.filter.clear()
            for f in record.filter: new_record.filter.add(f)
            
            if strand == '+':
                new_record.alleles = record.alleles
            else: # Reverse complement for negative strand
                new_record.ref = reverse_complement(record.ref)
                new_record.alts = tuple(reverse_complement(alt) for alt in record.alts) if record.alts else None
            
            for key, value in record.info.items(): new_record.info[key] = value
            new_record.info['OriginalGenomicChrom'] = original_chrom
            new_record.info['OriginalGenomicPos'] = genomic_pos

            for sample in record.samples:
                for key, value in record.samples[sample].items():
                    new_record.samples[sample][key] = value
            
            vcf_out.write(new_record)
            converted_count += 1

    vcf_in.close()
    vcf_out.close()
    print("\nConversion complete.")
    print(f"  - New records created: {converted_count} (one input may map to multiple transcripts)")
    print(f"  - Original records skipped: {skipped_count} (intronic/intergenic or no matching chromosome)")
    print(f"Output written to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert VCF from genomic to transcriptomic coordinates after removing a prefix from chromosome names."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file with genomic coordinates.")
    parser.add_argument("-g", "--gtf", required=True, help="Gene annotation GTF file.")
    parser.add_argument("-o", "--output", required=True, help="Path for the output VCF file.")
    parser.add_argument(
        "-p", "--prefix", default="CHROMOSOME_", help="Prefix to remove from CHROM names (default: 'CHROMOSOME_')."
    )
    
    args = parser.parse_args()
    
    # Need to import re inside the main execution block
    import re
    convert_vcf_coordinates(args.vcf, args.gtf, args.output, args.prefix)
