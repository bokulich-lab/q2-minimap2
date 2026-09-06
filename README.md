# q2-minimap2

![CI](https://github.com/bokulich-lab/q2-minimap2/actions/workflows/ci.yaml/badge.svg)
[![codecov](https://codecov.io/gh/bokulich-lab/q2-minimap2/graph/badge.svg?token=PSCAYJUP01)](https://codecov.io/gh/bokulich-lab/q2-minimap2)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

QIIME 2 plugin for sequence alignment using minimap2.

## Installation
We provide two options for installing q2-minimap2 via conda environment files, depending on your preferred QIIME2 distribution:


#### 1: Amplicon distribution
Use the environment file `q2-minimap2-qiime2-amplicon-2025.10.yml` to create a new conda environment.

```shell
conda env create -f q2-minimap2-qiime2-amplicon-2025.10.yml
conda activate q2-minimap2-amplicon
```

#### 2: Moshpit distribution
Use the environment file `q2-minimap2-qiime2-moshpit-2025.10.yml` to create a new conda environment based on the MOSHPIT distro. 

```shell
conda env create -f q2-minimap2-qiime2-moshpit-2025.10.yml
conda activate q2-minimap2-moshpit
```

## Provided Actions


1. **build-index**

    Build a Minimap2 index database from reference sequences.

2. **minimap2-search**

    Search for top hits in a reference database using alignment between the query sequences and reference database sequences using Minimap2. Returns a report of the top M hits for each query (where M=maxaccepts).

3. **filter-reads**

    This method aligns long-read sequencing data (from a FASTQ file) to a set of reference sequences, identifying sequences that match or do not match the reference within a specified identity percentage. The alignment is performed using Minimap2, and the results are processed using Samtools.

4. **extract-reads**

    This method aligns long-read sequencing data (from a FASTA file) to a set of reference sequences, identifying sequences that match or do not match the reference within a specified identity percentage. The alignment is performed using Minimap2, and the results are processed using Samtools.

5. **align**

    Align reads to a set of reference sequences and keep the alignment itself, rather than
    only the reads it selects. Produces a coordinate-sorted BAM per sample
    (`SampleData[AlignmentMap]`), which can be summarised with `alignment-stats` or handed
    to any other tool that reads BAM.

6. **alignment-stats**

    Summarize Minimap2 alignments per sample: read counts, mapping rate, mapping quality
    and reference coverage. The output is metadata, so it works directly with
    `qiime metadata tabulate` and can be joined to a sample metadata file.

7. **classify-consensus-minimap2**

    Assign taxonomy to query sequences using Minimap2. Performs alignment between query and reference reads, then assigns consensus taxonomy to each query sequence.


<br>



### Examples

* build-index
  - Build Minimap2 index database
  ```shell
  qiime minimap2 build-index --i-reference reference.qza --o-index index.qza --verbose
  ```

<br>

* minimap2-search
  - Generate both hits and no hits for each query. Keep a maximum of one hit per query (primary).
  ```shell
  qiime minimap2 minimap2-search --i-query fasta_reads.qza --i-index index.qza --o-search-results paf.qza --verbose
  ```

  - Generate only hits for each query. Keep a maximum of one hit per query (primary mappings).
  ```shell
  qiime minimap2 minimap2-search --i-query fasta_reads.qza --i-index index.qza --o-search-results paf_only_hits.qza --p-output-no-hits false --verbose
  ```

  - Generate only hits for each query, limiting the number of hits to a maximum of 3 per query. Ensure that each hit has a minimum similarity percentage of 90% to be considered valid.
  ```shell
  qiime minimap2 minimap2-search --i-query fasta_reads.qza --i-index index.qza --o-search-results paf_only_hits_ma3.qza --p-maxaccepts 3 --p-output-no-hits false --verbose
  ```

<br>

* filter-reads
  - Keep mapped (single-end reads)
  ```shell
  qiime minimap2 filter-reads --i-query single-end-reads.qza --i-index index.qza --o-filtered-query mapped_se.qza --verbose
  ```

  - Keep unmapped (single-end reads)
  ```shell
  qiime minimap2 filter-reads --i-query single-end-reads.qza --i-index index.qza --p-keep unmapped --o-filtered-query unmapped_se.qza --verbose
  ```

  - Keep mapped (paired-end reads)
  ```shell
  qiime minimap2 filter-reads --i-query paired-end-reads.qza --i-index index.qza --o-filtered-query mapped_pe.qza --verbose
  ```

  - Keep mapped reads with mapping percentage >= 98% (paired-end reads)
  ```shell
  qiime minimap2 filter-reads --i-query paired-end-reads.qza --i-index index.qza --p-min-per-identity 0.98  --o-filtered-query mapped_pe_over_98p_id.qza --verbose
  ```

<br>

* extract-reads
  - Extract mapped
  ```shell
  qiime minimap2 extract-reads --i-sequences fasta_reads.qza --i-index index.qza --o-extracted-reads mapped_fasta.qza --verbose
  ```
  - Extract unmapped
  ```shell
  qiime minimap2 extract-reads --i-sequences fasta_reads.qza --i-index index.qza --p-extract unmapped --o-extracted-reads unmapped_fasta.qza --verbose
  ```
  - Extract mapped reads with mapping percentage >= 87%
  ```shell
  qiime minimap2 extract-reads --i-sequences fasta_reads.qza --i-index index.qza --p-min-per-identity 0.87 --o-extracted-reads mapped_fasta_ido_ver_87.qza --verbose
  ```

<br>

* align
  - Align reads and keep the BAM
  ```shell
  qiime minimap2 align --i-query single-end-reads.qza --i-index index.qza --o-alignment alignment.qza --verbose
  ```
  - Use every available core
  ```shell
  qiime minimap2 align --i-query single-end-reads.qza --i-index index.qza --p-n-threads auto --o-alignment alignment.qza
  ```

<br>

* alignment-stats
  - Summarize an alignment, then view it as a table
  ```shell
  qiime minimap2 alignment-stats --i-alignments alignment.qza --o-stats alignment-stats.qza
  qiime metadata tabulate --m-input-file alignment-stats.qza --o-visualization alignment-stats.qzv
  ```

<br>

* filter-reads also reports how many reads it removed per sample
  ```shell
  qiime minimap2 filter-reads --i-query single-end-reads.qza --i-index index.qza \
    --o-filtered-query mapped_se.qza --o-filter-stats filter-stats.qza
  ```

<br>

* classify-consensus-minimap2
  - Assign taxonomy to query sequences using Minimap2
  ```shell
  qiime minimap2 classify-consensus-minimap2 --i-query n1K_initial_reads_SILVA132.fna.qza --i-index ccm_index.qza --i-reference-taxonomy raw_taxonomy.qza --p-n-threads 8 --output-dir classification_output --verbose
  ```
