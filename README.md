# TRAM-LOG
TRAM-LOG: a statistical framework for trans-ancestry association mapping by leveraging the local genetic structure

## Installation
``` shell
git clone https://github.com/YangLabHKUST/TRAM-LOGbeta.git
cd TRAM
conda env create -f environment.yml
conda activate tram
```
check the installation status
```shell
python ./src/TRAM.py -h
```

## Quick start

### Data preparation

Input files of TRAM include:

- GWAS summay statistics files of the target and auxiliary populations
- LDscore files (from LDXR)

The TRAM format GWAS summary statistics file has at least 11 fields:

- SNP: SNP rsid
- CHR: chromosome
- BP: base pair
- A1: effect allele
- A2: other allele
- FRQ: effect allele frequency
- BETA: marginal effect size
- SE: standard error
- N: sample size
- Z: Z-scores
- P: p-value 


### Usage
``` shell
python <install path>/src/TRAM.py \\
        --out test \\
        --sumstats-popu1 ss_file1,ss_name1 \\
        --sumstats-popu2 ss_file2,ss_name2 \\
        --ldscores ldsc_annot_1mb_TGP_hm3_chr@_std \\
        --out-harmonized 
```

## Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.



