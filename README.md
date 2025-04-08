# azmigueldario/bwa_align_stats

## Introduction

**azmigueldario/bwa_align_stats** is a bioinformatics pipeline that evaluates the taxonomic classification of isolates (in this case Giardia duodenalis) by aligning against reference genomes. Besides, it will use a high quality reference to produce a `pileup` format file and do variant calling through different approaches.

Its input are short read `fastq` files and reference genomes (`fasta`) as well as parameters to determine optional steps and flags for every module.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,organism
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz,
SAMPLE_PAIRED_END,/path/to/fastq/files/AEG588A1_S1_L002_R1_001.fastq.gz,/path/to/fastq/files/AEG588A1_S1_L002_R2_001.fastq.gz,"Pseudomonas_aeruginosa"
```

Each row represents a fastq file (paired end) with an accompanying sample_id and organism for mapping.

Now, you can run the pipeline using:

```bash
nextflow run azmigueldario/bwa_align_stats \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --ref_genome reference_genome_v01.fa.gz \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

azmigueldario/bwa_align_stats was originally written by azmigueldario.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Jimmy Liu](https://github.com/jimmyliu1326)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
