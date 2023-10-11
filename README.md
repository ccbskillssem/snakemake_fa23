# CCB Skills Seminar (Fall 2023)
**A gentle introduction to `snakemake`**: A tool for automating and streamlining your analyses

Created by [Stacy Li](stacy.li) for the [Center for Computational Biology](ccb.berkeley.edu) at UC Berkeley.

--------------

## What is this?
This is a repository containing a sample configuration + workflow for learning `snakemake`. The accompanying slides will be linked after the live workshop.

## Table of contents

- [Requirements](#requirements)
- [Setup](#setup)
- [Getting started](#getting-started)

## Requirements
You need `conda` (Anaconda distribution, `conda`, and `miniconda`) or a `conda`-like (`mamba`, `micromamba`, etc) package manager installed to run this tutorial. That's it! All other packages will be installed into their own isolated environments as you go along.

## Setup
Use`git clone` or download a copy of this repository to a suitable location.

(I will add a link for Datahub support here after I figure that out.)

## Getting started
Set this repository to your working directory, then run the following commands:

```
conda env create --file envs/snakemake_fa23.yml
conda activate snakemake_fa23
```

You will need to tweak this if you use a `conda`-like manager with a different command alias.