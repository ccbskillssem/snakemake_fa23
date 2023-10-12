# CCB Skills Seminar (Fall 2023)
## **A gentle introduction to `snakemake`**: A tool for automating and streamlining your analyses 🐍

Created by [Stacy Li](https://linktr.ee/stacy_li) for the [Center for Computational Biology](ccb.berkeley.edu) at UC Berkeley.

--------------

## What is this?
This is a repository containing a sample configuration + workflow for learning `snakemake`. I've tried to make this a minimally fussy example of how `snakemake` works, and all the good it does for me in my work 😊

## Table of contents

- [Slides](#slides)
- [Requirements](#requirements)
- [Live setup](#live-setup)
- [Independent setup](#independent-setup)
- [Commands](#commands)
- [Credits](#credits)
- [Contact info](#contact-info)

## Slides

Slides are included with the repo: click [here](https://github.com/ccbskillssem/snakemake_fa23/blob/main/ccb_snakemake_fa23_slidesv2.pdf) for a web view via GitHub. **Note that the slides have been updated with corrections made following the live workshop.**

## Requirements
You need `conda` (Anaconda distribution, `conda`, and `miniconda`) or a `conda`-like (`mamba`, `micromamba`, etc) package manager installed to run this tutorial. That's it! All other packages will be installed into their own isolated environments as you go along.

**Don't have this yet?** Click [here](https://docs.conda.io/projects/miniconda/en/latest/) and select the appropriate installer, or use one of the DataHub setup options below.

## Live setup
**This section refers to the setup commands used during the live (1 hour) workshop. If you are working through this on your own, please go to (#independent-setup).** 

The setup below pre-populates the `conda` environments you'll need for the workshop: this will save us precious time while we go over the lecture.
### Local machine or savio interactive (`srun`)
Use `git clone` or download a copy of this repository to a suitable location. This will be your working directory.
```
conda env create --file envs/snakemake_fa23.yml
conda activate snakemake_fa23
snakemake --cores all --use-conda --conda-create-envs-only output/visuals/vcf_heatmap.pdf
```

### DataHub
This is a UC Berkeley only option – you need a Calnet ID to proceed.
First, click [this link]() and proceed through the various pages until you launch your JupyterHub instance. If everything went well, you should have a cloned `snakemake_fa23` folder.
Open up a terminal instance, then run the commands below:
```
cd snakemake_fa23
conda env create --file envs/snakemake_fa23.yml
bash
conda activate snakemake_fa23
snakemake --cores all --use-conda --conda-create-envs-only output/visuals/vcf_heatmap.pdf
```

## Independent setup
This section refers to setup commands for independently working through the materials. `conda` environments will be generated on-the-fly as you run the workflow.

### Local machine or savio interactive (`srun`)
Use `git clone` or download a copy of this repository to a suitable location. This will be your working directory.
```
conda env create --file envs/snakemake_fa23.yml
conda activate snakemake_fa23
```

### DataHub
This is a UC Berkeley only option – you need a Calnet ID to proceed.
First, click [this link](https://biology.datahub.berkeley.edu/user/stacy-l/git-pull?repo=https://github.com/ccbskillssem/snakemake_fa23) and proceed through the various pages until you launch your JupyterHub instance. If everything went well, you should have a cloned `snakemake_fa23` folder.
Open up a terminal instance, then run the commands below:
```
cd snakemake_fa23
conda env create --file envs/snakemake_fa23.yml
bash
conda activate snakemake_fa23
```

## Commands
This is a list of commands from the live workshop. All of the rules are documented with comments as well. All examples below will run the workflow to generate the `output/visuals/vcf_heatmap.pdf` target file.

To do a dry run of the workflow, where `n` = max number of jobs to run in parallel:
```
snakemake -pj{n} --use-conda output/visuals/vcf_heatmap.pdf -np
```

For example, to run with 10 maximum jobs:

```
snakemake -pj10 --use-conda output/visuals/vcf_heatmap.pdf -np
```

To perform a *real* run of workflow:
```
snakemake -pj{n} --use-conda output/visuals/vcf_heatmap.pdf
```

To create a visualization of the rule graph:
```
snakemake --rulegraph --use-conda output/visuals/vcf_heatmap.pdf | dot -Tpng > output/visuals/workflow_rulegraph.png
```

To create a visualization of the *full* DAG (directed acyclic graph of every sample processed):
```
snakemake --dag --use-conda output/visuals/vcf_heatmap.pdf | dot -Tpng > output/visuals/workflow_dag.png
```

Once you feel comfortable with how the above commands work, I recommend trying out a few more things:
* Change the sample table from `config/subset.tsv` to `config/all_samples.tsv`. How can you evaluate what jobs might need to be run to update the vcf heatmap, without running the actual jobs?
* Actually run the workflow to update the vcf heatmap. Make sure it looks different!
* Create the full DAG for the workflow of all samples. It's probably quite visually busy.
* Try to implement a `rule all` in the `Snakefile` to ensure that both the coverage plot of your choice (`.html` or `.pdf`) and the vcf heatmap are always up to date.

## Credits and thanks
The short-read *Fructilactobacillus sanfranciscensis* data used in this workshop is from [Rogalski et al 2020](https://www.sciencedirect.com/science/article/pii/S0944501320304936). The `vcfR` heatmap plotting script is a modified version of a script from [Olawoye et al 2020](https://doi.org/10.7717/peerj.10121).

Thank you to my wonderful research group, [The Sudmant Lab](sudmantlab.org). As always, I am especially for my advisor Dr. Peter Sudmant and mentor Dr. Juan Manuel Vazquez. Everything I do is only possible because I stand upon the shoulders of giants. 🌟

## Contact info
If you'd like to stay in touch, please feel free to connect with me using any of the platforms [here](https://linktr.ee/stacy_li). If you're local, come and say hi at the next CCB event 👋