# LOG-TRAM
Leveraging the local genetic structure for trans-ancestry association mapping

## Installation
``` shell
$ git clone https://github.com/YangLabHKUST/LOG-TRAM.git
$ cd LOG-TRAM
$ conda env create -f environment.yml
$ conda activate tram
```
check the installation status
```shell
$ python ./src/LOG-TRAM.py -h
usage: LOG-TRAM.py [-h] --out OUT --sumstats-popu1 FILE,PHENOTYPE [FILE,PHENOTYPE ...] --sumstats-popu2 FILE,PHENOTYPE
                   [FILE,PHENOTYPE ...] --ldscores LDSCORES [--use_snps USE_SNPS] [--out-harmonized] [--out-reg-coef]
                   [--reg-int-ident] [--reg-int-diag] [--allowed-chr-values ALLOWED_CHR_VALUES [ALLOWED_CHR_VALUES ...]]
                   [--remove-palindromic-snps] [--num_threads NUM_THREADS]

Leverage local genetic architecture for trans-ancestry association mapping

optional arguments:
  -h, --help            show this help message and exit
  --out OUT             output file path
  --sumstats-popu1 FILE,PHENOTYPE [FILE,PHENOTYPE ...]
                        summary statisitcs F(file path),P(phenotype) of population 1, separated by whitespace
  --sumstats-popu2 FILE,PHENOTYPE [FILE,PHENOTYPE ...]
                        summary statisitcs F(file path),P(phenotype) of population 2, separated by whitespace
  --ldscores LDSCORES   specifies prefix of the LD score files computed by S-LDXR (popu1 <corresponding to population of --sumstats-
                        popu1>, popu2 <corresponding to population of --sumstats-popu2>, trans-ethnic), If the filename prefix
                        contains the symbol @, LOG-TRAM will replace the @ symbol with chromosome number, then add the suffix
                        _pop1.gz/_pop2.gz/_te.gz
  --use_snps USE_SNPS   SNPs list file (one rsID per line), If specified, this list will be used to restrict the final list of SNPs
                        reported
  --out-harmonized      If specified, LOG-TRAM will output harmonized summary statistics to disk
  --out-reg-coef        If specified, LOG-TRAM will output LD score regression coeficients to disk
  --reg-int-ident       Optional argument indicating that the LDscore regression intercept matrix should be set to be the identity
                        matrix.
  --reg-int-diag        Optional argument indicating that the LDscore regression intercept matrix should have off-diagonal elements
                        set to zero
  --allowed-chr-values ALLOWED_CHR_VALUES [ALLOWED_CHR_VALUES ...]
                        specify the allowed values for the chromosome
  --remove-palindromic-snps
                        This option removes the SNPs whose major and minor alleles form a base pair
  --num_threads NUM_THREADS
                        number of threads
```

## Reproducibility

We provide source codes and datasets for reproducing the experiments of LOG-TRAM meta-analysis of 29 EAS and EUR traits, and 17 AFR and EUR/EAS traits in the `demos` directory.
+ [BMI](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/BMI-demo.ipynb)
+ [SBP](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/SBP-demo.ipynb)
+ [T2D](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/T2D-demo.ipynb)
+ [26 Other Traits](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/26OtherTraits-demo.ipynb)
+ [17 AFR Traits](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/AFR-Traits-demo.ipynb)
+ [Simulations in Supplementary text](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/Simulation_simplified.ipynb) and [codes for generating figures](https://github.com/YangLabHKUST/LOG-TRAM/blob/main/demos/Plot_simplified.ipynb)

## Quick start

We illustrate the usage of LOG-TRAM by applying it to the GWAS summary statistics of BMI from BBJ male and UKBB with 1 Mbp non-overlapping sliding windows as local regions. The GWAS datasets and LDscores files involved in the following example are availabel from [here](https://www.dropbox.com/sh/9asugdlu1lbal8o/AAB0martsgaBoR8B4hq2pc25a?dl=0)

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
1       752566  rs3094315       G       A       0.8438  -0.0035920490000000004  0.006645739     -0.5405041      0.5889  85894
1       846808  rs4475691       C       T       0.1411  0.004290303     0.0069305669999999995   0.6190406999999999      0.5359  85894
1       854250  rs7537756       A       G       0.1766  0.0013868110000000002   0.0063270719999999996   0.2191868       0.8265  85894
1       861808  rs13302982      A       G       0.5404  0.0138975119697173      0.00484123214245955     2.8706559736788 0.004093        85894
1       863124  rs4040604       G       T       0.5421  0.01388188      0.0048426       2.866618        0.004145        85894
1       880238  rs3748592       A       G       0.9454  -0.01863968     0.01061939      -1.75525        0.07931 85894
1       882803  rs2340582       A       G       0.9462  -0.0187335737094156     0.0106935272024001      -1.75186104218362       0.07986 85894
1       884815  rs4246503       A       G       0.9422  -0.01772236     0.01033875      -1.7141680000000001     0.08656 85894
1       888659  rs3748597       T       C       0.9526  -0.01932025     0.01135427      -1.701585       0.0888  85894
```

LDscore files were computed by [S-LDXR](https://github.com/huwenboshi/s-ldxr) with easily accessible 1000 Genomes project genotypes as reference panels. 
For reproducibility, we provide the LDscore files of EUR, EAS, AFR, and trans-ancestries for 1 Mbp non-overlapping sliding windows [here](https://www.dropbox.com/sh/9asugdlu1lbal8o/AAB0martsgaBoR8B4hq2pc25a?dl=0)


### Usage
Once the input files are formatted, LOG-TRAM will automatically preprocess the datasets, including SNPs overlapping and minor allele matching. It takes 8 mins to run the following meta-analysis for the whole genome (computing environment: 20 CPU cores of Intel(R) Xeon(R) Gold 6230N CPU @ 2.30GHz processor, 1TB of memory, and a 22 TB solid-state disk). 

``` shell
python <install path>/src/LOG-TRAM.py \
        --out BMI_meta \
        --sumstats-popu1 BMI_harmonized_pop1_UKB.txt,BMI_UKB \
        --sumstats-popu2 BMI_harmonized_pop2_BBJ.txt,BMI_BBJ \
        --ldscores ./LDscoresEUR-EAS/ldsc_annot_EUR_EAS_1mb_TGP_hm3_chr@_std
```

### Results

LOG-TRAM will output two meta-analysis files, corresponding to EAS and EUR respectively. LOG-TRAM will add the inputed phenotype name after `--out` argument automatically. Usually, we focus on the under-represented populations such as EAS:

``` shell
$ head BMI_meta_TRAM_pop2_BMI_BBJ.txt
CHR     BP      SNP     A1      A2      FRQ     BETA    SE      Z       P       N       N_eff
1       752566  rs3094315       G       A       0.8438  -0.0015213296260925097  0.0028724631954416806   -0.5296254547340108     5.963716420681353e-1    85894   142023.18563988016
1       846808  rs4475691       C       T       0.1411  0.0024243137461232656   0.0032506101165059613   0.7458026829526787      4.557866190070761e-1    85894   142023.18563988016
1       854250  rs7537756       A       G       0.1766  0.0018498950124382176   0.0030984860402657115   0.597031901515225       5.504860818628083e-1    85894   142023.18563988016
1       861808  rs13302982      A       G       0.5404  0.009935609855423763    0.0035296626853253818   2.8148893368001393      4.879403304370678e-3    85894   142023.18563988016
1       863124  rs4040604       G       T       0.5421  0.010249240519351707    0.0035059696552997184   2.923368290954452       3.4626668356417944e-3   85894   142023.18563988016
1       880238  rs3748592       A       G       0.9454  -0.016950700361663937   0.006774302735757054    -2.502205912970559      1.234221174580461e-2    85894   142023.18563988016
1       882803  rs2340582       A       G       0.9462  -0.01697948824519616    0.006711745870831338    -2.5298169167857703     1.1412205878599073e-2   85894   142023.18563988016
1       884815  rs4246503       A       G       0.9422  -0.016311334186057062   0.006506282892858627    -2.507012752851644      1.2175631870410657e-2   85894   142023.18563988016
1       888659  rs3748597       T       C       0.9526  -0.017798097900681556   0.007087384241213514    -2.5112364865424786     1.2030907045241506e-2   85894   142023.18563988016
```
**N** is the original GWAS sample size, **N_eff** is the computed effective sample size. **N_eff** should be larger than **N** as LOG-TRAM can brorrow information from the large-scale auxiliary dataset.


## Contact information

Please contact Jiashun Xiao (jxiaoae@connect.ust.hk), Mingxuan Cai (mcaiad@ust.hk) or Prof. Can Yang (macyang@ust.hk) if any enquiry.



