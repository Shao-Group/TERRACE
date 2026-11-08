# TERRACE 2.0

TERRACE 2.0 is a circRNA assembler for paired-end RNA-seq data. It builds on its predecessor [TERRACE](https://github.com/Shao-Group/TERRACE) and reconstructs circRNA transcripts directly from read alignments.

A circRNA forms when the 3’ end of a downstream exon joins back to the 5’ end of an upstream exon via a back-splicing junction (BSJ). Reads that cross this BSJ appear in an alignment as either chimeric or outward reads. TERRACE 2.0 collects both kinds of evidence and assembles each candidate into a complete circRNA with its exons. Compared with the original TERRACE, it makes use of outward reads and filters candidates with learned models built into the program.

The input is a coordinate-sorted BAM file and the reference genome. The output is a GTF file with one `circRNA` entry per assembled transcript, followed by one `exon` entry for each of its exons.

## Installation

TERRACE 2.0 is built from source. Start by cloning the repository, then follow [INSTALL.md](INSTALL.md) to set up the dependencies (e.g., htslib, Boost, and zlib) and compile the program:

```
git clone -b outward https://github.com/Shao-Group/TERRACE.git
cd TERRACE
```

## Usage

```
terrace -i <input.bam> -o <output.gtf> -fa <reference-genome.fa> [--read_length <read-length>] [-r reference_annotation.gtf] [options]
```

| Argument | Required | Description |
| --- | --- | --- |
| `-i` | yes | Alignment file from an RNA-seq aligner such as STAR. It must be coordinate-sorted. |
| `-o` | yes | Output GTF file written by the original TERRACE (see "Where the output goes"). |
| `-fa` | yes | Reference genome in FASTA format (indexed with `samtools faidx`). |
| `--read_length` | optional | Length of the sequenced reads. Note that the default is 100 bps, so specify this value explicitly. |
| `-r` | optional | Reference gene annotation in GTF format. |
| `--library_type` | optional | One of `empty`, `unstranded`, `first`, `second`. These correspond to the usual Illumina `fr-unstranded`, `fr-firststrand` and `fr-secondstrand` libraries. If left empty, TERRACE infers it from the alignments. |
| `--results_tag` | optional | Name of the directory where results are written. Defaults to `results_YYYYMMDD_HHMMSS`. |
| `--help` | | Print the user manual and exit. |

A few notes on the inputs:

- Align your reads with STAR and sort the result. If the BAM file is not sorted yet, run `samtools sort input.bam -o input.sort.bam` first.
- Use the same reference genome for the FASTA that you aligned against. GRCh37/hg19 or GRCh38 from GENCODE or UCSC both work. Index the FASTA once with `samtools faidx reference-genome.fa`.

### Where the output goes

The GTF file you pass to `-o` is the output of the original TERRACE. It is written in the folder you run the program from and lists the assembled circRNAs. The flag is kept so that one run gives you both the original TERRACE's circRNAs and TERRACE 2.0's, which makes it easy to compare the two on the same input.

In addition, TERRACE 2.0 creates a results directory in that folder, named `results_<date>_<time>` unless you choose a different name with `--results_tag`. It contains `terrace_2.0.gtf`, the circRNAs assembled by the TERRACE 2.0 pipeline. This is the file we score in the evaluation section below.

The results directory also holds a few supporting files. You can ignore them unless you want to see how a particular candidate was handled:

| File | Contents |
| --- | --- |
| `circrna_paths.txt` | The selected circRNAs, grouped by bundle (a cluster of reads at one locus) and connected component. |
| `circrna_features.tsv` | One row per output circRNA, with the evidence features used to judge it (e.g., remaining chimeric and outward read support, minimum edge weight, minimum vertex weight). |
| `circrna_diagnostics.tsv` | A log of every candidate the assembler considered, including the ones it rejected. For each rejected candidate, a column names the filter that removed it. |
| `circrna_candidate_pool.tsv`, `circrna_selection_pool.tsv` | All candidate circRNAs generated from the reads, and the smaller set that the greedy path selection step chose from. |

## Running the small example

A small example is included in the `example` directory. `example-input.bam` contains aligned paired-end reads (150 bps) from chromosome 21 of the GRCh37/hg19 reference genome, so it is quick to run and is a good way to confirm that your build works.

Before running it, make sure you have built TERRACE 2.0 (see [INSTALL.md](INSTALL.md)) and have the matching reference genome. If you do not have the genome yet, you can download it (e.g., hg19) from UCSC:

```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
gunzip hg19.fa.gz
samtools faidx hg19.fa
```

Then run TERRACE 2.0 from the `example` directory:

```
cd example
../src/terrace -i example-input.bam -o example-output.gtf --read_length 150 -fa /path/to/hg19.fa
```

When it finishes, there will be a new `results_<date>_<time>` directory inside `example`. Open `terrace_2.0.gtf` there to see the reconstructed circRNAs. On this example, you should get 4 circRNAs.

## Evaluating accuracy

TERRACE 2.0 is evaluated by comparing its circRNAs against a ground-truth set and reporting recall and precision. A predicted circRNA is counted as correct only if it exactly matches a true one, meaning it is on the same chromosome and has identical coordinates for every exon.

For this, you need two GTF files: the `terrace_2.0.gtf` file from your TERRACE 2.0 run, and a ground-truth GTF listing the circRNAs you expect to find. The ground truth should use the same layout as TERRACE 2.0's output (with each `circRNA` line followed by its `exon` lines).

To compare TERRACE 2.0's results against the ground truth, run the script below. It only requires Python, with no extra packages to install:

```
python3 evaluate/compare_chains_gtf_two.py <ground_truth.gtf> <terrace_2.0.gtf>
```

The script prints a short summary. These are the lines to focus on particularly:

```
len of groups in groups1: ...
len of groups in groups2: ...
matching items: ...
Percentage of match against group 1/Recall: ...
Percentage of match against group 2/Precision: ...
```

The first two lines give the number of circRNAs in the ground truth and in TERRACE 2.0, respectively; `matching items` is how many of them appear in both (the true positives). Recall and precision are then computed and reported as percentages in the following lines.