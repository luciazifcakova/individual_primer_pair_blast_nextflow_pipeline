# Individual Primer-Pair BLAST Nextflow Pipeline

This repository contains a **sanitized excerpt of a production-oriented Nextflow DSL2 workflow** that I developed for long-read amplicon analysis.

It is published as an example of my Nextflow workflow design and is **not intended to be a directly runnable distribution of the production pipeline**. Company-specific modules, helper scripts, databases, configuration data, and other proprietary components are intentionally not included.

## What the workflow does

The complete workflow processes long-read amplicon sequencing data through multiple configurable stages, including:

- optional sample demultiplexing
- read QC and filtering
- primer-aware preparation of sequences for analysis
- parallelized BLAST analysis
- taxonomic classification
- aggregation of QC and analysis results
- generation of customer-facing and internal reports

Different stages can be enabled or disabled depending on the analysis required for a particular project.

## Nextflow design

The workflow is implemented in **Nextflow DSL2** and combines nf-core components with custom modules.

The orchestration shown in `main.nf` demonstrates several aspects of my workflow development experience:

- modular workflow composition
- Nextflow channels and metadata propagation
- handling of sample, project, and assay metadata
- conditional execution of analysis stages
- splitting computationally expensive analyses into independently scheduled tasks
- regrouping and joining outputs from parallel processes
- execution on HPC/Slurm
- use of Nextflow caching and `-resume` for efficient recovery and reproducibility
- separation of workflow logic from computational resource configuration

One important design goal was to make large BLAST analyses scalable. BLAST workloads are divided into independent Nextflow tasks that can be executed in parallel by Slurm and cached independently. This allows failed or modified analyses to be resumed without unnecessarily repeating completed computation.

## Validation and reproducibility

For production bioinformatics workflows, successful execution alone is not sufficient. I therefore treat validation as part of workflow development.

This includes checking:

- read retention after QC/filtering
- expected intermediate files and sample counts
- taxonomic assignment behaviour
- effects of filtering and classification thresholds
- handling of empty or unusual samples
- consistency of aggregated results
- computational resource usage and parallelization behaviour

I use small controlled datasets and inspection of intermediate outputs when developing or modifying workflow components.

## AI-assisted development

I use AI tools extensively to accelerate software development, particularly for drafting implementations, debugging, investigating unfamiliar software behaviour, comparing possible approaches, and improving documentation.

I do not treat generated code as authoritative.

I define the biological and computational requirements, review the proposed implementation, test it on controlled inputs, inspect logs and intermediate outputs, and determine whether the resulting behaviour is scientifically and computationally correct.

AI therefore helps me iterate substantially faster, while responsibility for **workflow architecture, biological assumptions, validation criteria, and acceptance of the final implementation remains with me**.

## Open-source contribution

I have also contributed to the **nf-core/modules** project, including work on the BLAST `makeblastdb` module:

https://github.com/nf-core/modules/blob/master/modules/nf-core/blast/makeblastdb/meta.yml

This provided experience working with nf-core conventions and contributing to a community-maintained bioinformatics codebase.

## Repository contents

- `main.nf` — main Nextflow DSL2 workflow and orchestration logic
- `nextflow.config` — Nextflow configuration and computational resource settings

The custom modules and supporting components referenced by the workflow are part of the original development environment and are intentionally not included in this public repository.
