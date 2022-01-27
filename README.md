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

We illustrate the usage of LOG-TRAM using the GWAS summary statistics of BMI from BBJ and UKBB. The datasets and LDscores files involved in the following example is availabel from [here](https://www.dropbox.com/sh/9asugdlu1lbal8o/AAB0martsgaBoR8B4hq2pc25a?dl=0)

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
$ head BMI_harmonized_pop1_UKB.txt
CHR     BP      SNP     A1      A2      FRQ     BETA    SE      Z       P       N
1       846808  rs4475691       T       C       0.197757        0.00279009838327684     0.00659822032318736     0.42285620161423204     0.67    72390
1       854250  rs7537756       G       A       0.208596        0.00330351434039922     0.006468351398041789    0.510719677567196       0.61    72390
1       880238  rs3748592       A       G       0.053943899999999996    -0.00685935129704364    0.0116336527874558      -0.589612860411297      0.56   72390
1       882803  rs2340582       A       G       0.0537767       -0.00751214563103344    0.0116506946664748      -0.644780920458747      0.52    72390
1       884815  rs4246503       A       G       0.053765099999999996    -0.007067130889541141   0.0116518800153349      -0.606522799774817      0.54   72390
1       888659  rs3748597       T       C       0.0539094       -0.00693882584725204    0.0116371625581102      -0.596264408321446      0.55    72390
1       891945  rs13303106      A       G       0.345632        -0.00238214650323602    0.00552621420712329     -0.431063005151959      0.67    72390
1       918573  rs2341354       A       G       0.409366        -0.0010370878670859     0.00534478952201353     -0.19403717635923       0.85    72390
1       928836  rs9777703       C       T       0.0347591       -0.00176196795009865    0.0143480735553437      -0.122801708766152      0.9     72390
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
1       846808  rs4475691       T       C       0.197757        -0.0005827564488972137  0.003308918373456196    -0.17611690078909958    8.602021035447315e-1    72390   59086.33115808354
1       854250  rs7537756       G       A       0.208596        0.0004537332696092035   0.0034152970863215158   0.13285323593851744     8.943094512798043e-1    72390   59086.33115808354
1       880238  rs3748592       A       G       0.0539438999999999      -0.006995947452271588   0.006259761974268365    -1.1176059858233296     2.6373531115192668e-1   72390   59086.33115808354
1       882803  rs2340582       A       G       0.0537767       -0.006944825916337652   0.006063426216986064    -1.1453633090945246     2.520586906286916e-1    72390   59086.33115808354
1       884815  rs4246503       A       G       0.0537650999999999      -0.006733842734770233   0.006082743865309179    -1.1070403232288588     2.6827649490206893e-1   72390   59086.33115808354
1       888659  rs3748597       T       C       0.0539094       -0.006682554860854716   0.006060609842431025    -1.102620864004375      2.701918458104707e-1    72390   59086.33115808354
1       891945  rs13303106      A       G       0.345632        -0.0029196306422206433  0.0024276676481119796   -1.2026484121462293     2.29112404945011e-1     72390   59086.33115808354
1       918573  rs2341354       A       G       0.409366        -0.0006359488661905613  0.0028145389695380697   -0.22595134516646434    8.212392626195399e-1    72390   59086.33115808354
1       928836  rs9777703       C       T       0.0347591       -0.0025919938637620405  0.008687291813913964    -0.2983661559072512     7.654237170043164e-1    72390   59086.33115808354
```
**N_eff** is the computed effective sample size 


## Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.



