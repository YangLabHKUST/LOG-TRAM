# LOG-TRAM
Leveraging the local genetic structure for trans-ancestry association mapping

## Installation
``` shell
$ git clone https://github.com/YangLabHKUST/LOG-TRAM.git
$ cd TRAM
$ conda env create -f environment.yml
$ conda activate tram
```
check the installation status
```shell
$ python ./src/LOG-TRAM.py -h
usage: LOG-TRAM.py [-h] --out OUT --sumstats-popu1 FILE,PHENOTYPE
               [FILE,PHENOTYPE ...] --sumstats-popu2 FILE,PHENOTYPE
               [FILE,PHENOTYPE ...] --ldscores LDSCORES
               [--annot-names ANNOT_NAMES] [--use_snps USE_SNPS]
               [--out-harmonized] [--out-reg-coef {yes,no}]
               [--local-anno {yes,no}] [--reg-int-ident] [--reg-int-diag]
               [--reg-ld-coef REG_LD_COEF]
               [--allowed-chr-values ALLOWED_CHR_VALUES [ALLOWED_CHR_VALUES ...]]
               [--remove-palindromic-snps] [--num_threads NUM_THREADS]

Leverage local genetic architecture for trans-ancestry association mapping

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             output file path
  --sumstats-popu1 FILE,PHENOTYPE [FILE,PHENOTYPE ...]
                        summary statisitcs F(file path),P(phenotype)
                        of population 1, separated by whitespace
  --sumstats-popu2 FILE,PHENOTYPE [FILE,PHENOTYPE ...]
                        summary statisitcs F(file path),P(phenotype)
                        of population 2, separated by whitespace
  --ldscores LDSCORES   specifies prefix of the LD score files computed by
                        S-LDXR (popu1 <corresponding to population of
                        --sumstats-popu1>, popu2 <corresponding to population
                        of --sumstats-popu2>, trans-ethnic), If the filename
                        prefix contains the symbol @, LOG-TRAM will replace the @
                        symbol with chromosome number, then add the suffix
                        _pop1.gz/_pop2.gz/_te.gz
  --use_snps USE_SNPS   SNPs list file (one rsID per line), If specified, this
                        list will be used to restrict the final list of SNPs
                        reported
  --out-harmonized      If specified, LOG-TRAM will output harmonized summary
                        statistics to disk
  --out-reg-coef {yes,no}
                        If specified, LOG-TRAM will output LD score regression
                        coeficients to disk
  --local-anno {yes,no}
                        Optional argument indicating that the annotation is
                        non-overlapping sliding windows across the genome. If
                        this arugument is specified, then the LDScore
                        regression will be perfromed on the 'base' LD score
                        column and one local genetic annotation at a time to
                        avoid singular matrix problem (because the 'base' LD
                        score column almost equals to the sum of other LD
                        score columns in this case).
  --annot-names ANNOT_NAMES
                        Functinoal annotation name list file, one name per
                        line (do not include "base")
  --reg-int-ident       Optional argument indicating that the LDscore
                        regression intercept matrix should be set to be the
                        identity matrix.
  --reg-int-diag        Optional argument indicating that the LDscore
                        regression intercept matrix should have off-diagonal
                        elements set to zero
  --reg-ld-coef REG_LD_COEF
                        Optional argument indicating the file (dict data type
                        store in .npy file) containing the regression
                        coeficient for LDscores of functional annotations
  --allowed-chr-values ALLOWED_CHR_VALUES [ALLOWED_CHR_VALUES ...]
                        specify the allowed values for the chromosome
  --remove-palindromic-snps
                        This option removes the SNPs whose major and minor
                        alleles form a base pair
  --num_threads NUM_THREADS
                        number of threads
```

## Quick start

We illustrate the usage of LOG-TRAM by applying it to the GWAS summary statistics of BMI from BBJ and UKBB with 1 Mbp non-overlapping sliding windows as local regions. The GWAS datasets and LDscores files involved in the following example are availabel from [here](https://www.dropbox.com/sh/9asugdlu1lbal8o/AAB0martsgaBoR8B4hq2pc25a?dl=0)

### Data preparation

Input files of LOG-TRAM include:

- GWAS summay statistics files of the target and auxiliary populations
- LDscore files (from [S-LDXR](https://github.com/huwenboshi/s-ldxr))

The LOG-TRAM format GWAS summary statistics file has at least 11 fields:

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
e.g.,
``` shell
$ head BMI_harmonized_pop2_BBJ.txt
CHR     BP      SNP     A1      A2      FRQ     BETA    SE      Z       P       N
1       846808  rs4475691       T       C       0.1411  -0.001356749    0.002191696     -0.6190406999999999     0.5359  85894
1       854250  rs7537756       G       A       0.1766  -0.0004385598   0.00200085      -0.2191868      0.8265  85894
1       880238  rs3748592       A       G       0.9454  -0.005894544    0.003358235     -1.75525        0.07931 85894
1       882803  rs2340582       A       G       0.9462  -0.00592423480012795    0.00338168077117785     -1.75186104218362       0.07986 85894
1       884815  rs4246503       A       G       0.9422  -0.005604451999999999   0.003269488     -1.714168       0.08656 85894
1       888659  rs3748597       T       C       0.9526  -0.006109762    0.003590631     -1.701585       0.0888  85894
1       891945  rs13303106      A       G       0.6436  -0.00305593142187254    0.0015930828132597302   -1.91825019794141       0.0551  85894
1       918573  rs2341354       A       G       0.8325  -0.0005534582070046848  0.0020432226342925104   -0.270875135051706      0.7865  85894
1       928836  rs9777703       C       T       0.9876  -0.003776945    0.006894676     -0.547806       0.5839  85894
```

LDscore files were computed by [S-LDXR](https://github.com/huwenboshi/s-ldxr) with easily accessible 1000 Genomes project genotypes as reference panels. 
For reproducibility, we provide the LDscore files of EUR, EAS, and trans-ancestry for 1 Mbp non-overlapping sliding windows [here](https://www.dropbox.com/sh/9asugdlu1lbal8o/AAB0martsgaBoR8B4hq2pc25a?dl=0)


### Usage
It takes 8 mins to run the following meta-analysis for the whole genome (computing environment: 20 CPU cores of Intel(R) Xeon(R) Gold 6230N CPU @ 2.30GHz processor, 1TB of memory, and a 22 TB solid-state disk).

``` shell
python <install path>/src/LOG-TRAM.py \\
        --out BMI_meta \\
        --sumstats-popu1 BMI_harmonized_pop1_UKB.txt,BMI_UKB \\
        --sumstats-popu2 BMI_harmonized_pop2_BBJ.txt,BMI_BBJ \\
        --ldscores ldsc_annot_1mb_TGP_hm3_chr@_std 
```

### Results

LOG-TRAM will output two meta-analysis files, corresponding to EAS and EUR respectively. Usually, we focus on the under-represented populations such as EAS:

``` shell
$ head BMI_meta_TRAM_pop2_BMI_BBJ.txt
CHR     BP      SNP     A1      A2      FRQ     BETA    SE      Z       P       N       N_eff
1       846808  rs4475691       T       C       0.1411  -0.0011173325152100634  0.002191792209263744    -0.5097803115129204     6.1020538060194545e-1  85894    97954.28694930511
1       854250  rs7537756       G       A       0.1766  -0.00022805612123296662 0.002003734026852899    -0.11381556542769089    9.093839993491452e-1   85894    97954.28694930511
1       880238  rs3748592       A       G       0.9454  -0.006041632961525742   0.0033925199267538005   -1.7808688208080172     7.493388325577553e-2   85894    97954.28694930511
1       882803  rs2340582       A       G       0.9462  -0.006100310957504757   0.0034156438901755794   -1.7859915007682998     7.410063429224731e-2   85894    97954.28694930511
1       884815  rs4246503       A       G       0.9422  -0.005759375254588021   0.0033025343506386017   -1.7439259196423942     8.117200827667153e-2   85894    97954.28694930511
1       888659  rs3748597       T       C       0.9526  -0.006275907311551924   0.003625821062505641    -1.7308927283948559     8.347089881296858e-2   85894    97954.28694930511
1       891945  rs13303106      A       G       0.6436  -0.0030794980721019377  0.0016050267743185124   -1.9186583808917952     5.502758190554889e-2   85894    97954.28694930511
1       918573  rs2341354       A       G       0.8325  -0.0006012663641836504  0.0020557010902238363   -0.29248725266676834    7.699141048816635e-1   85894    97954.28694930511
1       928836  rs9777703       C       T       0.9876  -0.0038126379088130754  0.006975857051358909    -0.5465475970540948     5.8468957920787386e-1  85894    97954.28694930511
```
**N** is the original GWAS sample size, **N_eff** is the computed effective sample size 


## Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.



