import sys
import os
from Bio import SeqIO
import pandas as pd
import numpy as np
from collections import Counter


genomes_dirpath = sys.argv[1]
queries_dirpath = sys.argv[2]
output_dirpath = sys.argv[3]
blasted_dirpath = f'{output_dirpath}/blasted'


def filepaths_fasta(dirpath):
    files = os.listdir(dirpath)
    files.sort()
    fastas = list(filter(lambda file: ".fasta" in file, files))
    filepaths = [f'{dirpath}/{fasta}' for fasta in fastas]
    return filepaths

def extract_filename(filepath):
    file = os.path.basename(filepath)
    filename = ".".join(file.split(".")[:-1])
    return filename

def extract_scaffolds_ids(filepath):
    fasta_seqs = SeqIO.parse(open(genome_filepath), 'fasta')
    scaffolds = []
    for seq in fasta_seqs:
        scaffolds.append(seq.id)
    scaffolds.sort()
    return scaffolds

def read_blast_output(blasted_filepath):
    alignments = pd.read_table(blasted_filepath, \
                               skiprows = 5, \
                               header = None, \
                               names = ["query", "subject_accession", \
                                        "pident", "alen", \
                                        "mismatches", "gap_opens", \
                                        "q_start", "q_end", \
                                        "s_start", "s_end", \
                                        "evalue", "bit_score"], \
                               skipfooter = 1, \
                               engine = 'python')
    return alignments

def add_info_to_bed(bed):
    bed["s_start"] - 1
    bed["pident"] * 10
    bed['strand'] = np.where(bed['s_start'] > bed['s_end'], '-', '+')
    neg_strand_num = bed[bed['strand'] == '-'].shape[0]
    if neg_strand_num != 0:
        bed.loc[bed['strand'] == '-', ['s_start', 's_end']] = \
        bed.loc[bed['strand'] == '-', ['s_end', 's_start']].values
    return bed

def filter_bed(bed, e_val_thresh, len_thresh, query_len):
    bed_filtered = bed[bed['evalue'] < e_val_thresh]
    query_filter_len = int(query_len * len_thresh)
    bed_filtered = bed[bed['alen'] > query_filter_len]
    return bed_filtered

def write_tsv_line(bed, all_scaffolds, tsv_dirpath, genome_filename, guery_filename, e_val_thresh, len_thresh, query_len):
    scaffolds = bed["subject_accession"].values.tolist()
    query_per_scaf = Counter(scaffolds)
    query_per_scaf = sorted(query_per_scaf.items())

    query_per_scaf_all = []
    for scaffold in all_scaffolds:
        query_per_scaf_all.append([scaffold, 0])

    j = 0
    for i in range(len(query_per_scaf_all)):
        if j >= len(query_per_scaf):
            break
        if query_per_scaf_all[i][0] == query_per_scaf[j][0]:
            query_per_scaf_all[i][1] = query_per_scaf[j][1]
            j += 1

    query_num = []
    for pair in query_per_scaf_all:
        query_num.append(str(pair[1]))
    
    all_scaffolds_string = "\t".join(all_scaffolds)
    query_num_string = "\t".join(query_num)
    
    tsv_filepath = f'{tsv_dirpath}/{genome_filename}.tsv'
    with open(tsv_filepath, 'a+') as tsv:
        if os.stat(tsv_filepath).st_size == 0:
            all_scaffolds_string = "\t".join(all_scaffolds)
            tsv.write(f'query\tquery_len\te_value_threshold\talign_len_threshold\t{all_scaffolds_string}\n')
        tsv.write(f'{query_filename}\t{query_len}\t{e_val_thresh}\t{len_thresh}\t{query_num_string}\n')


genomes_filepaths = filepaths_fasta(genomes_dirpath)
queries_filepaths = filepaths_fasta(queries_dirpath)


e_val_threshs = [0.01, 1e-15]
len_threshs = [0.75, 0.90, 0.95]


for genome_filepath in genomes_filepaths:
    genome_filename = extract_filename(genome_filepath)
    all_scaffolds = extract_scaffolds_ids(genome_filepath)
    
    for query_filepath in queries_filepaths:
        query_filename = extract_filename(query_filepath)
        blasted_filename = f'{genome_filename}_{query_filename}'
        blasted_file = f'{blasted_filename}.blasted'
        blasted_filepath = f'{blasted_dirpath}/{blasted_file}'
        
        print()
        print()
        print(f'### Converting {blasted_file} to .bed format ###')
        bed = read_blast_output(blasted_filepath)
        bed = add_info_to_bed(bed)
        
        filter_params = '10_0'
        
        bed_file = f'{blasted_filename}.bed'
        bed_dirpath = f'{output_dirpath}/{filter_params}/bed'
        os.makedirs(bed_dirpath, exist_ok = True)
        bed_filepath = f'{bed_dirpath}/{bed_file}'
        print()
        print(f'Writing unfiltered .bed file')
        bed[["subject_accession", "s_start", "s_end", "query", "pident", "strand"]].to_csv(bed_filepath, sep="\t", header=False, index=False)

        query = SeqIO.parse(query_filepath, 'fasta')
        query_len = len(next(query).seq)
        
        tsv_dirpath = f'{output_dirpath}/{filter_params}/tsv'
        os.makedirs(tsv_dirpath, exist_ok = True)
        write_tsv_line(bed, all_scaffolds, tsv_dirpath, genome_filename, query_filename, 10, 0, query_len)

        for e_val_thresh in e_val_threshs:
            for len_thresh in len_threshs:
                print()
                print(f'Filtering {bed_file} with {e_val_thresh} e-value and {len_thresh} alignment length')
                bed_filtered = filter_bed(bed, e_val_thresh, len_thresh, query_len)
                bed_filtered = bed_filtered[["subject_accession", "s_start", "s_end", "query", "pident", "strand"]]
                
                filter_params = f'{e_val_thresh}_{int(len_thresh * 100)}'
                
                bed_filtered_dirpath = f'{output_dirpath}/{filter_params}/bed'
                os.makedirs(bed_filtered_dirpath, exist_ok = True)
                bed_filtered_file = f'{blasted_filename}.bed'
                bed_filtered_filepath = f'{bed_filtered_dirpath}/{bed_filtered_file}'
                print()
                print(f'Writing {bed_filtered_file} to {filter_params} directory')
                bed_filtered.to_csv(bed_filtered_filepath, sep="\t", header=False, index=False)
                
                tsv_dirpath = f'{output_dirpath}/{filter_params}/tsv'
                os.makedirs(tsv_dirpath, exist_ok = True)
                write_tsv_line(bed_filtered, all_scaffolds, tsv_dirpath, genome_filename, query_filename, e_val_thresh, len_thresh, query_len)
                print()
    
    print(f'All .tsv files written for {genome_filename}')


print('\n\n\n\n')
print('##### Drawing heatmaps #####')


dirs = [dir_ for dir_ in os.listdir(output_dirpath) if os.path.isdir(f'{output_dirpath}/{dir_}')]
dirs.remove('blasted')
dirs.sort()

for dir_ in dirs:
    print()
    print()
    print(f'Processing {dir_} directory')
    heatmap_dirpath = f'{output_dirpath}/{dir_}'
    tsv_dirpath = f'{heatmap_dirpath}/tsv'
    
    command = f'python ./1.1.1_draw_heatmap.py {tsv_dirpath} {heatmap_dirpath}'
    os.system(command)
    