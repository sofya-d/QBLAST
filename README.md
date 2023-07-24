# QBLAST: Quantitative analysis of genomic features with NCBI BLAST
### What does this pipeline do?
This pipeline aligns input query sequences of genomic features on input genomes and produces tables and plot with query sequences number in genomes.
### Output data
* BLAST output files with default filtering parameters
* [BED](https://genome.ucsc.edu/FAQ/FAQformat.html) files filtered with all combinations of alignment length and BLAST e-value parameter values (0.01, 1e-15 and 0.75, 0.90, 0.95 respectively)
* Summary [TSV](https://en.wikipedia.org/wiki/Tab-separated_values) tables for each combination of parameters where rows are query sequences, columns are scaffolds/contigs, and values in cells are number of specific query sequence in specific scaffold/contig.
* Heat map plot in [SVG](https://en.wikipedia.org/wiki/SVG) format were rows are input genomes, columns are input query sequences, and values in cells are number of specific query sequence in whoe genome.
### Prerequisites
* Unix OS<br>
  This pipeline was tested with Ubuntu 20.04.3 LTS.
* [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) command line tool<br>
  This pipeline was tested with version 2.13.0+. NCBI blast must be added to [PATH variable](https://en.wikipedia.org/wiki/PATH_(variable)).
* [Python](https://www.python.org/)<br>
  This pipeline was tested with version 3.10.10.
* [Biopython](https://biopython.org/) Python module<br>
  This pipeline was tested with version 1.79.
* [Matplotlib](https://matplotlib.org/) Python module<br>
  This pipeline was tested with version 3.7.1.
* [NumPy](https://numpy.org/) Python module<br>
  This pipeline was tested with version 1.24.3.
* [pandas](https://pandas.pydata.org/) Python module<br>
  This pipeline was tested with version 1.5.3.
* [seaborn](https://seaborn.pydata.org/) Python module<br>
  This pipeline was tested with version 0.12.2.
### Installation
Download and inflate archive with pipeline via GitHub GUI or if you have Git installed paste the following command in shell:
```shell
git clone https://github.com/sofya-d/QBLAST.git
```
### Input
1. Path to directory with genome [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) files. Files must have '.fasta' extension
2. Path to directory where BLAST databases will be written
3. Path to directory with query sequence FASTA files. Files must have '.fasta' extension
4. Path to directory with script files of this pipeline
5. Path to output directory
### Pipeline running example
```shell
./1_run_blast.sh ./genomes ./databases ./queries ./ ./output
```
