# SnakeAlign
SGE snakemake pipeline implementation for genome alignment and coverage statistics analysis.  
Currently only works with pair end sequencing data.  
** This pipeline was thought for S. cerevisiae thus, the 2m sequence analysis.

## **Download:**

```
git clone --recursive https://github.com/Iseez/SnakeAlign.git
```

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

## **Dependencies:**
- [Bwa](https://bio-bwa.sourceforge.net)
- [Samtools](http://www.htslib.org)
- [Fastp](http://www.htslib.org)
- [Picard](https://broadinstitute.github.io/picard/)
- [R](https://www.r-project.org)
  - [ggplot2](https://ggplot2.tidyverse.org)
  - [data.table](https://rdatatable.gitlab.io/data.table/)
  - [zoo](https://cran.r-project.org/web/packages/zoo/index.html)

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
**Edit `cluster.yaml`:**
```YAML
__default__:
  M: 'your@email.com'
  N: 'rule_{rule}.{wildcards}'
  l: 'vf=4G'
  q: 'all.q@compute-00-12.cm'
  pe: 'openmp 4'
  mem: '6G'
  m: 'a'
  r: 'n'
  e: 'st/log-{rule}.{wildcards}.err'
  o: 'st/log-{rule}.{wildcards}.out'
  S: '/bin/bash'
```
#### **Ready to run**
To use Snakemake the module `python37/3.7.0` should be loaded on your current session of the cluster. Also you should add `r/3.6.1` to the list of modules that autoload (`module initadd r/3.6.1`).  
\**Other R modules should also work but weren't tested.

**Run `Snakefile`**  
To run the pipeline you can run the following line:
```
snakemake data/done.txt

```
But to use the full potential of the cluster the following line is recommended:
```
snakemake -p --latency-wait 60 -j 10 --cluster-config cluster.yaml --cluster 'qsub -V -S {cluster.S} -N {cluster.N} -cwd -e {cluster.e} -o {cluster.o} -q {cluster.q} -M {cluster.M} -l {cluster.l} -pe {cluster.pe} -m {cluster.m} -r {cluster.r}' data/done.txt
```
\**In the previous line the -j **10** denotes the maximum amount of samples to be processed at any given time,you should probably leave it in less than 20.
