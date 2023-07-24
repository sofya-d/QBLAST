#!/bin/bash

genomes_dirpath=$1
blast_dbs_dirpath=$2
queries_dirpath=$3
scripts_dirpath=$4
output_dirpath=$5

echo
echo
echo "##### Input #####"
echo "Path to genomes is $genomes_dirpath"
echo "Path to BLAST databases is $blast_dbs_dirpath"
echo "Path to queary sequences is $queries_dirpath"
echo "Path to script that is running now is $scripts_dirpath"
echo "Path to BLAST output is $output_dirpath"


mkdir -p $blast_dbs_dirpath
export BLASTDB=$blast_dbs_dirpath
mkdir -p $output_dirpath
echo


echo
echo "##### Creating BLAST databases #####"

cd $blast_dbs_dirpath
for genome_filepath in $genomes_dirpath/*.fasta
do
    genome_file=$(basename -- "$genome_filepath")
    genome_filename=${genome_file%.*}

    if [ -f ${genome_filename}.nhr ]; then
        echo "BLAST databases exist for genome $genome_filename"
    else 
        echo "Creating BLAST database for $genome_filename"
        
        makeblastdb -in ${genome_filepath} \
                -input_type 'fasta' \
                -dbtype nucl \
                -title $genome_filename \
                -out $genome_filename
    
    fi
done
echo


echo
echo "##### Blasting query sequences over genomes #####"
cd $scripts_dirpath
blasted_dirpath="$output_dirpath/blasted"
mkdir -p $blasted_dirpath

for query_filepath in $queries_dirpath/*.fasta
do
	for genome_filepath in $genomes_dirpath/*.fasta
    do
        
	    genome_file=$(basename -- "$genome_filepath")
	    genome_filename=${genome_file%.*}
        
        query_file=$(basename -- "$query_filepath")
        query_filename="${query_file%.*}"
        
        blasted_filepath=${blasted_dirpath}/${genome_filename}_${query_filename}.blasted
        
        echo
        echo "Blasting $query_filename over $genome_filename"
        
        blastn -query ${query_filepath} \
               -db ${genome_filename} \
               -task megablast \
               -out ${blasted_filepath} \
               -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
               -num_threads 32
        
    done
done
echo
echo


echo
echo
echo "##### Creating .bed and .tsv files out of BLAST output files #####"
python ./1.1_make_tables.py $genomes_dirpath $queries_dirpath $output_dirpath
