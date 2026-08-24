import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os
import sys

vcf_file = snakemake.input.vcf
genome1 = snakemake.input.genome1
genome2 = snakemake.input.genome2
insertion_fasta = snakemake.output.insertion_fasta
blast_tsv1 = snakemake.output.blast_results_genome1
blast_tsv2 = snakemake.output.blast_results_genome2
filtered_vcf = snakemake.output.filtered_vcf
stats_file = snakemake.output.blast_stats
min_length = snakemake.params.min_ins_length
min_pident = snakemake.params.min_pident
min_qcovs = snakemake.params.min_qcovs
threads = snakemake.threads

os.makedirs(os.path.dirname(insertion_fasta), exist_ok=True)
os.makedirs(os.path.dirname(filtered_vcf), exist_ok=True)


insertions = []
vcf_data = []

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt = fields[4]
        if len(alt) > len(ref):
            ins_seq = alt[len(ref):]
            if len(ins_seq) >= min_length:
                ins_id = f"{chrom}_{pos}_{len(ins_seq)}"
                record = SeqRecord(Seq(ins_seq), id=ins_id)
                insertions.append(record)
                vcf_data.append({'ins_id': ins_id, 'chrom': chrom, 'pos': pos, 'vcf_line': line.strip()})

SeqIO.write(insertions, insertion_fasta, "fasta")
vcf_df = pd.DataFrame(vcf_data)
print(f"Extracted {len(insertions)} insertions")

columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qcovs', 'qcovhsp']

cmd1 = ['blastn', '-query', insertion_fasta, '-subject', genome1, '-out', blast_tsv1, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp', '-max_target_seqs', '1', '-max_hsps', '1', '-num_threads', str(threads)]
print(f"Running: {' '.join(cmd1)}")
subprocess.run(cmd1, check=True)

cmd2 = ['blastn', '-query', insertion_fasta, '-subject', genome2, '-out', blast_tsv2, '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp', '-max_target_seqs', '1', '-max_hsps', '1', '-num_threads', str(threads)]
print(f"Running: {' '.join(cmd2)}")
subprocess.run(cmd2, check=True)


blast1_df = pd.DataFrame(columns=columns)
blast2_df = pd.DataFrame(columns=columns)

if os.path.getsize(blast_tsv1) > 0:
    blast1_df = pd.read_csv(blast_tsv1, sep='\t', names=columns)

if os.path.getsize(blast_tsv2) > 0:
    blast2_df = pd.read_csv(blast_tsv2, sep='\t', names=columns)

print(f"BLAST1 results: {len(blast1_df)} hits")
print(f"BLAST2 results: {len(blast2_df)} hits")


perfect_hits = set()

if not blast1_df.empty:
    perfect1 = blast1_df[(blast1_df['pident'] >= min_pident) & (blast1_df['qcovs'] >= min_qcovs)]
    perfect_hits.update(perfect1['qseqid'].tolist())

if not blast2_df.empty:
    perfect2 = blast2_df[(blast2_df['pident'] >= min_pident) & (blast2_df['qcovs'] >= min_qcovs)]
    perfect_hits.update(perfect2['qseqid'].tolist())

print(f"Found {len(perfect_hits)} insertions with perfect hits")
print(f"Perfect hits: {perfect_hits}")


removed_count = 0
kept_count = 0
filtered_lines = []

with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            filtered_lines.append(line)
        else:
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            if len(alt) > len(ref):
                ins_seq = alt[len(ref):]
                if len(ins_seq) >= min_length:
                    ins_id = f"{chrom}_{pos}_{len(ins_seq)}"
                    if ins_id not in perfect_hits:
                        filtered_lines.append(line)
                        kept_count += 1
                    else:
                        removed_count += 1
                else:
                    filtered_lines.append(line)
            else:
                filtered_lines.append(line)

with open(filtered_vcf, 'w') as f:
    f.writelines(filtered_lines)

print(f"Removed {removed_count} insertions with perfect hits")
print(f"Kept {kept_count} insertions")



with open(stats_file, 'w') as f:
    f.write("BLAST Filter Statistics\n")
    f.write("========================\n\n")
    f.write(f"Total insertions extracted: {len(insertions)}\n")
    f.write(f"BLAST1 hits: {len(blast1_df)}\n")
    f.write(f"BLAST2 hits: {len(blast2_df)}\n")
    f.write(f"Perfect hits (100% ID, 100% coverage): {len(perfect_hits)}\n")
    f.write(f"Insertions removed: {removed_count}\n")
    f.write(f"Insertions kept: {kept_count}\n")
    f.write(f"\nPerfect hit IDs:\n")
    for hit in sorted(perfect_hits):
        f.write(f"  {hit}\n")

print("Done!")
