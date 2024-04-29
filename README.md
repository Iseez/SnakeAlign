# SnakeAlign
SGE snakemake pipeline implementation for genome alignment and coverage statistics analysis.  
Currently only works with pair end sequencing data.  
** This pipeline was thought for S. cerevisiae thus, the 2m sequence analysis.

## **Download:**
---
```
git clone --recursive https://github.com/Iseez/SnakeAlign.git
```
---

## **Content:**
```
.
├── config.yaml
├── data
│   └── reads
├── LICENSE
├── README.md
├── scripts
│   ├── make_depth_plot.R
│   └── normalized_coverage.R
├── Snakefile
└── st
```
- `config.yaml`: Configuration file with necessary parameters to run the pipeline. ***Edit this file before running!***
- `data/reads`: put the pair end sequencing FASTQs here
  - FASTQs shouldn't be compressed and should be named as  *Sample1*\_R1.fatsq and *Sample1*\_R2.fatsq
- `scripts`: scripts
- `Snakefile`: Snakemake worflow
---
## **Dependencies:**
- [Bwa](https://bio-bwa.sourceforge.net)
- [Samtools](http://www.htslib.org)
- [Fastp](http://www.htslib.org)
- [Picard](https://broadinstitute.github.io/picard/)
- [R](https://www.r-project.org)
  - ggplot2
  - data.table
  - zoo
---
## **Usage:**
#### **Before starting**

Put the FASTQs corresponding to the samples to analyse inside `./data/reads/`.

**Edit `config.yaml`:**
```YAML
# Pipeline parameters

#Fasta with the reference.
reference_file: "/path/to/ref.fasta"

#Reference name.
reference_name: "Ref"

#Sequencing platform for picard
platform: "illumina|solid|etc"

#Quality for depth assesment
quality: 30

#Parameters for plot making. SIZE is the window size, STEP is the step size
size: 1000
step: 1000

#Samples to align
samples:
  - Sample1
  - Sample2

```  


**Run `Snakefile`**
To run the 
```
snakemake -p --latency-wait 60 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -q {cluster.q} -cwd -e {cluster.e} -o {cluster.o}' data/done.txt
```

---
