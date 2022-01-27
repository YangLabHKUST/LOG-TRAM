import pandas as pd
import numpy as np
import sys
import argparse
import os
from utils import *

__version__ = '1.0.1'
SOFTWARE_CORRESPONDENCE_EMAIL1 = 'jxiaoae@connect.ust.hk'
SOFTWARE_CORRESPONDENCE_EMAIL2 = 'mcaiad@connect.ust.hk' 

HEADER = """
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
<> 
<> TRAM: Leverage local genetic architecture and functional annotation for trans-ancestry association mapping
<> Version: %s
<> MIT License
<>
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
<> Software-related correspondence: %s or %s
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
<> example:
    python <install path>/src/TRAM.py \\
        --out test \\
        --sumstats-popu1 ss_file1,ss_name1 \\
        --sumstats-popu2 ss_file2,ss_name2 \\
        --ldscores ldsc_annot_1mb_TGP_hm3_chr@_std \\
        --out-harmonized 
        
<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>       
""" % (__version__, SOFTWARE_CORRESPONDENCE_EMAIL1, SOFTWARE_CORRESPONDENCE_EMAIL2)


if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description='Leverage local genetic architecture for trans-ancestry association mapping')
    # Input / output 
    parser.add_argument('--out', type=str, help='output file path', required=True)
    parser.add_argument('--sumstats-popu1', type=str, nargs="+", metavar="FILE,PHENOTYPE", 
        help='summary statisitcs triples F(file path),P(phenotype) of population 1, separated by whitespace',required=True)
    parser.add_argument('--sumstats-popu2', type=str, nargs="+", metavar="FILE,PHENOTYPE", 
        help='summary statisitcs triples F(file path),P(phenotype) of population 2, separated by whitespace',required=True)
    parser.add_argument('--ldscores', type=str,
        help='specifies prefix of the LD score files computed by S-LDXR \
            (popu1 <corresponding to population of --sumstats-popu1>, \
            popu2 <corresponding to population of --sumstats-popu2>, trans-ethnic), \
            If the filename prefix contains the symbol @, TRAM will replace the @ symbol \
            with chromosome number, then add the suffix _pop1.gz/_pop2.gz/_te.gz',required=True)
    parser.add_argument('--annot-names', type=str,
       help='Functinoal annotation name list file, one name per line (do not include "base")',required=False)
    parser.add_argument('--use_snps', type=str, 
        help='SNPs list file (one rsID per line), If specified, this list will be used to restrict the final list of SNPs reported')
    parser.add_argument('--out-harmonized', action="store_true", 
        help='If specified, TRAM will output harmonized summary statistics to disk')
    parser.add_argument('--out-reg-coef', type=str, default='yes', choices=['yes','no'],
        help='If specified, TRAM will output LD score regression coeficients to disk')
    # Regression option
    parser.add_argument("--local-anno", type=str, default='yes', choices=['yes','no'],
        help="Optional argument indicating that the annotation is non-overlapping sliding windows across the genome.\
            If this arugument is specified, then the LDScore regression will be perfromed on the 'base' LD score column and one local genetic annotation at a time to \
                avoid singular matrix problem (because the 'base' LD score column almost equals to the sum of other LD score columns in this case).")
    parser.add_argument("--reg-int-ident", action="store_true",
        help="Optional argument indicating that the LDscore regression intercept matrix should be set to be the identity matrix.")
    parser.add_argument("--reg-int-diag", action="store_true",
        help="Optional argument indicating that the LDscore regression intercept matrix should have off-diagonal elements set to zero")
    parser.add_argument('--reg-ld-coef', type=str,
       help='Optional argument indicating the file (dict data type store in .npy file) containing the regression \
           coeficient for LDscores of functional annotations')
    # Summary Statistics Filtering Options
    parser.add_argument("--allowed-chr-values", type=str.upper, nargs="+",
        help='specify the allowed values for the chromosome')
    parser.add_argument("--remove-palindromic-snps", action="store_true",
        help='This option removes the SNPs whose major and minor alleles form a base pair')
    parser.add_argument('--num_threads', type=str, help='number of threads', default="22")
    args = parser.parse_args()

    os.environ["OMP_NUM_THREADS"] = args.num_threads
    os.environ["OPENBLAS_NUM_THREADS"] = args.num_threads
    os.environ["MKL_NUM_THREADS"] = args.num_threads
    os.environ["VECLIB_MAXIMUM_THREADS"] = args.num_threads
    os.environ["NUMEXPR_NUM_THREADS"] = args.num_threads

    logger = configure_logging(args.out)

    logger.info(HEADER)
    logger.info("See full log at: %s\n", os.path.abspath(args.out+'.log'))
    logger.info("\nProgram executed via:\n%s\n", ' '.join(sys.argv).replace("--", " \\ \n\t--"))

    logger.info("Reading in the base LD Scores")
    df_ldsc_pop1 = pd.concat(
        (pd.read_csv(args.ldscores.replace('@',str(c))+'_pop1.gz',
        compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']] for c in range(1,23)),
        ignore_index=True)
    df_ldsc_pop1 = df_ldsc_pop1.dropna()
    df_ldsc_pop2 = pd.concat(
        (pd.read_csv(args.ldscores.replace('@',str(c))+'_pop2.gz',
        compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']] for c in range(1,23)),
        ignore_index=True)
    df_ldsc_pop2 = df_ldsc_pop2.dropna()
    df_ldsc_te = pd.concat(
        (pd.read_csv(args.ldscores.replace('@',str(c))+'_te.gz',
        compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']] for c in range(1,23)),
        ignore_index=True)
    df_ldsc_te = df_ldsc_te.dropna()
    logger.info("LD score matrix shape: {}x{}".format(df_ldsc_pop1.shape[0],df_ldsc_pop1.shape[1])) 

    # Check / read in summary stats and then QC
    logger.info("Reading in summary statistics")
    sumstats_pop1 = {}
    sumstats_pop1_names = []
    for ssfile in args.sumstats_popu1:
        sumstats_pop1_names.append(ssfile.split(',')[1])
        sumstats_pop1[ssfile.split(',')[1]] = pd.read_csv(ssfile.split(',')[0],
            sep='\t',compression='infer',dtype={'CHR':int}).drop_duplicates('SNP').dropna()
        logger.info("Sumstats {} shape: {}x{}".format(ssfile.split(',')[1],
            sumstats_pop1[ssfile.split(',')[1]].shape[0],sumstats_pop1[ssfile.split(',')[1]].shape[1])) 
    sumstats_pop2 = {}
    sumstats_pop2_names = []
    for ssfile in args.sumstats_popu2:
        sumstats_pop2_names.append(ssfile.split(',')[1])
        sumstats_pop2[ssfile.split(',')[1]] = pd.read_csv(ssfile.split(',')[0],
            sep='\t',compression='infer',dtype={'CHR':int}).drop_duplicates('SNP').dropna()
        logger.info("Sumstats {} shape: {}x{}".format(ssfile.split(',')[1],
            sumstats_pop2[ssfile.split(',')[1]].shape[0],sumstats_pop2[ssfile.split(',')[1]].shape[1]))
    P1n, P2n = len(sumstats_pop1_names), len(sumstats_pop2_names)

    # matching SNPs
    logger.info('Matching SNPs from LD scores and GWAS Sumstats')
    snps_lst = [set(df_ldsc_pop1['SNP'].values),set(df_ldsc_pop2['SNP'].values),set(df_ldsc_te['SNP'].values)]+\
        [set(sumstats_pop1[_]['SNP'].values) for _ in sumstats_pop1.keys()]+\
            [set(sumstats_pop2[_]['SNP'].values) for _ in sumstats_pop2.keys()]
    snps_common = snps_lst[0].intersection(*snps_lst)
    if args.use_snps is not None:
        logger.info('loading user-supplied SNPs from {}'.format(args.use_snps))
        snps_common_inp = np.loadtxt(args.use_snps,dtype=str)
        snps_common = list(snps_common & set(snps_common_inp))           
    logger.info('Number of SNPS in initial intersection of all sources {}'.format(len(snps_common)))
    df_ldsc_pop1 = df_ldsc_pop1.loc[df_ldsc_pop1['SNP'].isin(snps_common)].reset_index(drop=True)
    df_ldsc_pop2 = df_ldsc_pop2.loc[df_ldsc_pop2['SNP'].isin(snps_common)].reset_index(drop=True)
    df_ldsc_te = df_ldsc_te.loc[df_ldsc_te['SNP'].isin(snps_common)].reset_index(drop=True)
    for k,v in sumstats_pop1.items():
        sumstats_pop1[k] = v.loc[v['SNP'].isin(snps_common)].sort_values(['CHR','BP']).reset_index(drop=True)
    for k,v in sumstats_pop2.items():
        sumstats_pop2[k] = v.loc[v['SNP'].isin(snps_common)].sort_values(['CHR','BP']).reset_index(drop=True)

    # Removing ambiguous SNPs
    logger.info('Removing ambiguous SNPs...')
    for k,v in sumstats_pop1.items():
        v['A1_int'] = v['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
        v['A2_int'] = v['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    for k,v in sumstats_pop2.items():
        v['A1_int'] = v['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
        v['A2_int'] = v['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
    sumstat_ref = sumstats_pop1[sumstats_pop1_names[0]]
    sumstat_ref_allelesum = sumstat_ref['A1_int'].values+sumstat_ref['A2_int'].values
    if args.remove_palindromic_snps:
        snps_rm_idx = (np.array([(sumstat_ref_allelesum % 2 == 0)]+\
            [(v['A1_int'].values+v['A2_int'].values != sumstat_ref_allelesum) for v in sumstats_pop1.values()]+\
            [(v['A1_int'].values+v['A2_int'].values != sumstat_ref_allelesum) for v in sumstats_pop2.values()])).sum(axis=0)>0
    else:
        snps_rm_idx = (np.array([(v['A1_int'].values+v['A2_int'].values != sumstat_ref_allelesum) for v in sumstats_pop1.values()]+\
            [(v['A1_int'].values+v['A2_int'].values != sumstat_ref_allelesum) for v in sumstats_pop2.values()])).sum(axis=0)>0
    logger.info('Number of ambiguous SNPS have been removed: {}'.format(snps_rm_idx.sum()))
    df_ldsc_pop1 = df_ldsc_pop1.loc[~snps_rm_idx,:].reset_index(drop=True)
    df_ldsc_pop2 = df_ldsc_pop2.loc[~snps_rm_idx,:].reset_index(drop=True)
    df_ldsc_te = df_ldsc_te.loc[~snps_rm_idx,:].reset_index(drop=True)
    for k,v in sumstats_pop1.items():
        sumstats_pop1[k] = v.loc[~snps_rm_idx,:].reset_index(drop=True)
    for k,v in sumstats_pop2.items():
        sumstats_pop2[k] = v.loc[~snps_rm_idx,:].reset_index(drop=True)

    logger.info('Aligning the minor allele according to a designated reference summary statistics DataFrame: {}'.format(sumstats_pop1_names[0]))
    allele_ref = sumstats_pop1[sumstats_pop1_names[0]][['SNP','A1','A2']].copy()
    allele_ref.columns = ['SNP','A1_ref','A2_ref']
    for k,v in sumstats_pop1.items():
        v = v.merge(allele_ref,left_on='SNP',right_on='SNP',how='inner')
        index_flip = v['A1'].values!=v['A1_ref'].values
        logger.info('{} SNPs need to be flipped for sumstats: {}'.format(index_flip.sum(),k))
        v.loc[index_flip,'BETA'] = -v.loc[index_flip,'BETA']
        v.loc[index_flip,'Z'] = -v.loc[index_flip,'Z']
        v['CHR'] = v['CHR'].astype(str)
        sumstats_pop1[k] = v
        if args.out_harmonized:
            v = v[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA','SE','Z','P','N']]
            v.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N']
            v.to_csv(args.out+'_harmonized_pop1_'+k+'.txt',sep='\t',index=None)
    for k,v in sumstats_pop2.items():
        v = v.merge(allele_ref,left_on='SNP',right_on='SNP',how='inner')
        index_flip = v['A1'].values!=v['A1_ref'].values
        logger.info('{} SNPs need to be flipped for sumstats: {}'.format(index_flip.sum(),k))
        v.loc[index_flip,'BETA'] = -v.loc[index_flip,'BETA']
        v.loc[index_flip,'Z'] = -v.loc[index_flip,'Z']
        v['CHR'] = v['CHR'].astype(str)
        sumstats_pop2[k] = v
        if args.out_harmonized:
            v = v[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA','SE','Z','P','N']]
            v.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N']
            v.to_csv(args.out+'_harmonized_pop2_'+k+'.txt',sep='\t',index=None)

    # convert to stdized units 
    for k,v in sumstats_pop1.items():
        conversion_factor_col = np.sqrt(2.0 * v['FRQ'] * (1.0 - v['FRQ']))
        v['BETA'] = v['BETA'] * conversion_factor_col
        v['SE'] = v['SE'] * conversion_factor_col
    for k,v in sumstats_pop2.items():
        conversion_factor_col = np.sqrt(2.0 * v['FRQ'] * (1.0 - v['FRQ']))
        v['BETA'] = v['BETA'] * conversion_factor_col
        v['SE'] = v['SE'] * conversion_factor_col

    # Run LD score regressions
    logger.info("Running LD Score regression for the whole genome.")
    Omega_global, Sigma_LD = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1,sumstats_pop2_names,\
        sumstats_pop2,df_ldsc_pop1[['base']].values,df_ldsc_pop2[['base']].values,df_ldsc_te[['base']].values)
    Omega_global = Omega_global.squeeze()
    logger.info("Regression coefficients (LD):\n%s", Omega_global)
    logger.info("Regression coefficients (Intercept):\n%s", Sigma_LD)

    if args.reg_int_diag:
        Sigma_LD = np.diag(np.diag(Sigma_LD))
        logger.info("Modified regression coefficients (Intercept):\n%s", Sigma_LD)
    if args.reg_int_ident:
        Sigma_LD = np.eye(Sigma_LD.shape[0])
        logger.info("Modified regression coefficients (Intercept):\n%s", Sigma_LD)

    # beta, se, N_eff
    tram_res_beta = np.zeros((df_ldsc_pop1.shape[0],P1n+P2n))
    tram_res_beta_se = np.zeros((df_ldsc_pop1.shape[0],P1n+P2n))

    if args.allowed_chr_values is None:
        chrs = [str(_) for _ in range(1,23)]
    else:
        chrs = args.allowed_chr_values

    if args.local_anno == 'yes':
        res_omegas_local = {}

        # Run LD score regressions for annotation
        logger.info("\nRunning LD Score regression and TRAM for local regions chromosome by chromosome.\n")

        for c in chrs:
            ldscore1_c = pd.read_csv(args.ldscores.replace('@',c)+'_pop1.gz',compression='infer',sep='\t')
            ldscore2_c = pd.read_csv(args.ldscores.replace('@',c)+'_pop2.gz',compression='infer',sep='\t')
            ldscorete_c = pd.read_csv(args.ldscores.replace('@',c)+'_te.gz',compression='infer',sep='\t')
            ldscore1_c = ldscore1_c.loc[ldscore1_c['SNP'].isin(df_ldsc_pop1['SNP'])]
            ldscore2_c = ldscore2_c.loc[ldscore2_c['SNP'].isin(df_ldsc_pop1['SNP'])]
            ldscorete_c = ldscorete_c.loc[ldscorete_c['SNP'].isin(df_ldsc_pop1['SNP'])]
            chr_snp_idx = df_ldsc_pop1.loc[df_ldsc_pop1['CHR']==c].index.values
            snp_num_c = chr_snp_idx.shape[0]
            logger.info("Processing chromosome {}, SNPs number: {}".format(c,chr_snp_idx.shape[0]))
            sumstats_pop1_c = {}
            beta_array = np.zeros((snp_num_c,P1n+P2n))
            se_array = np.zeros((snp_num_c,P1n+P2n))
            for i,ss_name in enumerate(sumstats_pop1_names):
                v = sumstats_pop1[ss_name].copy()
                v = v.loc[v['CHR']==c].reset_index(drop=True)
                sumstats_pop1_c[ss_name] = v
                beta_array[:,i] = v['BETA'].values
                se_array[:,i] = v['SE'].values
            sumstats_pop2_c = {}
            for i,ss_name in enumerate(sumstats_pop2_names):
                v = sumstats_pop2[ss_name].copy()
                v = v.loc[v['CHR']==c].reset_index(drop=True)
                sumstats_pop2_c[ss_name] = v
                beta_array[:,i+P1n] = v['BETA'].values
                se_array[:,i+P1n] = v['SE'].values

            name_windows = list(ldscore1_c.columns.values[4:])
            n_window = len(name_windows)
            ldscore1_ck = ldscore1_c[['base']+name_windows].values
            ldscore2_ck = ldscore2_c[['base']+name_windows].values
            ldscorete_ck = ldscorete_c[['base']+name_windows].values
            logger.info("find {} local regions in chromosome {}".format(len(name_windows),c))
            omegas_local = np.zeros((1+len(name_windows),P1n+P2n,P1n+P2n))
            omegas_local_baseLD = np.zeros((len(name_windows),P1n+P2n,P1n+P2n))
            for i,name_window in enumerate(name_windows):
                if ldscore1_ck[:,i+1].max()<1 or ldscore2_ck[:,i+1].max()<1: # in case of emmpty window
                    omegas_local_k = np.zeros((2,P1n+P2n,P1n+P2n))
                    omegas_local_k[0] = Omega_global
                else:
                    omegas_local_k, _ = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1_c,sumstats_pop2_names,\
                        sumstats_pop2_c,ldscore1_ck[:,[0,i+1]],ldscore2_ck[:,[0,i+1]],ldscorete_ck[:,[0,i+1]],Sigma_LD)
                    if np.all(np.diag(omegas_local_k[1]) > 1e-10):
                        omegas_local_k[1] = make_positive_definite(omegas_local_k[1])
                omegas_local[i+1] = omegas_local_k.sum(axis=0)
                omegas_local_baseLD[i] = omegas_local_k[0]
            omegas_local_baseLD_mean = omegas_local_baseLD.mean(axis=0)
            omegas_local[0] = 2*omegas_local_baseLD_mean
            omegas_local -= omegas_local_baseLD_mean
            logger.info("Local regions regression coefficients (LD):\n%s \n...", omegas_local[:2])
            res_omegas_local[c] = dict(zip(['base']+name_windows,omegas_local))

            logger.info("Creating omega and sigma matrices for each SNP.")
            omega_ldscore_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))
            sigma_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))
            for j in range(snp_num_c):
                sigma_snps[j] = np.diag(se_array[j])@Sigma_LD@np.diag(se_array[j])
                ld_matrx_j = np.zeros_like(omegas_local)
                ld_matrx_j[:,:P1n,:P1n] = np.repeat(ldscore1_ck[j],P1n**2).reshape(ldscore1_ck[j].shape[0],P1n,P1n)
                ld_matrx_j[:,:P1n,P1n:] = np.repeat(ldscorete_ck[j],P1n*P2n).reshape(ldscorete_ck[j].shape[0],P1n,P2n)
                ld_matrx_j[:,P1n:,:P1n] = np.repeat(ldscorete_ck[j],P1n*P2n).reshape(ldscorete_ck[j].shape[0],P2n,P1n)
                ld_matrx_j[:,P1n:,P1n:] = np.repeat(ldscore2_ck[j],P2n**2).reshape(ldscore2_ck[j].shape[0],P2n,P2n)
                omega_ldscore_snps[j] = (omegas_local*ld_matrx_j).sum(axis=0)
            logger.info("Checking omega and sigma for validity based on positive (semi-)definiteness")
            fails_snps_omega = []
            fails_snps_sigma = []
            for j in range(snp_num_c):
                if not np.all(linalg.eigvalsh(omega_ldscore_snps[j]) >= 0.0):
                    fails_snps_omega.append(j)
                    omega_ldscore_snps[j] = make_positive_definite(omega_ldscore_snps[j])
                if not np.all(linalg.eigvalsh(sigma_snps[j]) >= 0.0):
                    fails_snps_sigma.append(j)
                    sigma_snps[j] = make_positive_definite(sigma_snps[j])
            logger.info("Adjusted %s SNPs to make omega positive definite.",len(fails_snps_omega))
            logger.info("Adjusted %s SNPs to make sigma positive definite.",len(fails_snps_sigma))

            # Run the TRAM method
            logger.info("Run the core TRAM method\n")
            new_betas, new_beta_ses = run_tram_method(beta_array,omega_ldscore_snps,sigma_snps)
            tram_res_beta[chr_snp_idx] = new_betas
            tram_res_beta_se[chr_snp_idx] = new_beta_ses

        if args.out_reg_coef == 'yes':
            logger.info("Saving local regions regression coefficients to file %s \n",args.out+'_TRAM_reg_coefs.npy')
            np.save(args.out+'_TRAM_reg_coefs.npy',res_omegas_local)
    
    else:
        logger.info("Running LD Score regression for functional annotations.")

        if args.annot_names is None:
            annot_names = pd.read_csv(args.ldscores.replace('@',str(22))+'_pop1.gz',
            compression='infer',sep='\t').columns.values[4:].tolist()
        else:
            annot_names = np.loadtxt(args.annot_names,dtype=str).tolist()
            if type(annot_names) == str:
                annot_names = [annot_names]
        
        n_window = len(annot_names)
        logger.info("Detecting {} functional annotations: {}...".format(n_window,','.join(annot_names[:3])))
        
        df_ldsc_pop1_all = pd.concat(
            (pd.read_csv(args.ldscores.replace('@',str(c))+'_pop1.gz',
            compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']+annot_names] for c in range(1,23)),
            ignore_index=True)
        df_ldsc_pop1_all = df_ldsc_pop1_all.loc[df_ldsc_pop1_all['SNP'].isin(df_ldsc_pop1['SNP'])]

        df_ldsc_pop2_all = pd.concat(
            (pd.read_csv(args.ldscores.replace('@',str(c))+'_pop2.gz',
            compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']+annot_names] for c in range(1,23)),
            ignore_index=True)
        df_ldsc_pop2_all = df_ldsc_pop2_all.loc[df_ldsc_pop2_all['SNP'].isin(df_ldsc_pop1['SNP'])]

        df_ldsc_te_all = pd.concat(
            (pd.read_csv(args.ldscores.replace('@',str(c))+'_te.gz',
            compression='infer',sep='\t',dtype={'CHR':str})[['CHR','SNP','BP','base']+annot_names] for c in range(1,23)),
            ignore_index=True)
        df_ldsc_te_all = df_ldsc_te_all.loc[df_ldsc_te_all['SNP'].isin(df_ldsc_pop1['SNP'])]

        ldscore1_all = df_ldsc_pop1_all[['base']+annot_names].values
        ldscore2_all = df_ldsc_pop2_all[['base']+annot_names].values
        ldscorete_all = df_ldsc_te_all[['base']+annot_names].values

        omegas_local, _ = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1,sumstats_pop2_names,\
            sumstats_pop2,ldscore1_all,ldscore2_all,ldscorete_all,Sigma_LD)

        for i,_ in enumerate(['base']+annot_names):
            if np.all(np.diag(omegas_local[i]) > 1e-10):
                omegas_local[i] = make_positive_definite(omegas_local[i])
        
        if args.reg_ld_coef is not None:
            inp_omegas = np.load(args.reg_ld_coef,allow_pickle=True).item()
            logger.info("Loading %s regression coeficent(s) from %s \n", len(inp_omegas.keys()), args.reg_ld_coef)
            for i,an in enumerate(['base']+annot_names):
                if an in inp_omegas.keys():
                    omegas_local[i] = inp_omegas[an]

        logger.info("Functional annotations regression coefficients (LD):\n%s \n...", omegas_local[:4])

        if args.out_reg_coef == 'yes':
            logger.info("Saving functional annotations regression coefficients to file %s \n",args.out+'_TRAM_reg_coefs.npy')
            np.save(args.out+'_TRAM_reg_coefs.npy',dict(zip(['base']+annot_names,omegas_local)))

        logger.info("\nRun the TRAM method chromosome by chromosome.\n")
        for c in chrs:
            chr_snp_idx = df_ldsc_pop1.loc[df_ldsc_pop1['CHR']==c].index.values
            snp_num_c = chr_snp_idx.shape[0]
            logger.info("Processing chromosome {}, SNPs number: {}".format(c,chr_snp_idx.shape[0]))
            beta_array = np.zeros((snp_num_c,P1n+P2n))
            se_array = np.zeros((snp_num_c,P1n+P2n))
            N_array = np.zeros((snp_num_c,P1n+P2n))
            for i,ss_name in enumerate(sumstats_pop1_names):
                v = sumstats_pop1[ss_name].copy()
                v = v.loc[v['CHR']==c]
                beta_array[:,i] = v['BETA'].values
                se_array[:,i] = v['SE'].values
                N_array[:,i] = v['N'].values
            for i,ss_name in enumerate(sumstats_pop2_names):
                v = sumstats_pop2[ss_name].copy()
                v = v.loc[v['CHR']==c]
                beta_array[:,i+P1n] = v['BETA'].values
                se_array[:,i+P1n] = v['SE'].values
                N_array[:,i+P1n] = v['N'].values
            
            logger.info("Creating omega and sigma matrices for each SNP.")
            omega_ldscore_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))
            sigma_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))

            for j in range(snp_num_c):
                sigma_snps[j] = np.diag(se_array[j])@Sigma_LD@np.diag(se_array[j])
                ld_matrx_j = np.zeros_like(omegas_local)
                ld_matrx_j[:,:P1n,:P1n] = np.repeat(ldscore1_all[chr_snp_idx[j]],P1n**2).reshape(ldscore1_all[chr_snp_idx[j]].shape[0],P1n,P1n)
                ld_matrx_j[:,:P1n,P1n:] = np.repeat(ldscorete_all[chr_snp_idx[j]],P1n*P2n).reshape(ldscorete_all[chr_snp_idx[j]].shape[0],P1n,P2n)
                ld_matrx_j[:,P1n:,:P1n] = np.repeat(ldscorete_all[chr_snp_idx[j]],P1n*P2n).reshape(ldscorete_all[chr_snp_idx[j]].shape[0],P2n,P1n)
                ld_matrx_j[:,P1n:,P1n:] = np.repeat(ldscore2_all[chr_snp_idx[j]],P2n**2).reshape(ldscore2_all[chr_snp_idx[j]].shape[0],P2n,P2n)
                omega_ldscore_snps[j] = (omegas_local*ld_matrx_j).sum(axis=0)
            logger.info("Checking omega and sigma for validity based on positive (semi-)definiteness")
            fails_snps_omega = []
            fails_snps_sigma = []
            for j in range(snp_num_c):
                if not np.all(linalg.eigvalsh(omega_ldscore_snps[j]) >= 0.0):
                    fails_snps_omega.append(j)
                    omega_ldscore_snps[j] = make_positive_definite(omega_ldscore_snps[j])
                if not np.all(linalg.eigvalsh(sigma_snps[j]) >= 0.0):
                    fails_snps_sigma.append(j)
                    sigma_snps[j] = make_positive_definite(sigma_snps[j])
            logger.info("Adjusted %s SNPs to make omega positive definite.",len(fails_snps_omega))
            logger.info("Adjusted %s SNPs to make sigma positive definite.",len(fails_snps_sigma))

            # Run the TRAM method
            logger.info("Run the core TRAM method\n")
            new_betas, new_beta_ses = run_tram_method(beta_array,omega_ldscore_snps,sigma_snps)
            tram_res_beta[chr_snp_idx] = new_betas
            tram_res_beta_se[chr_snp_idx] = new_beta_ses


    logger.info("Preparing results")
    sumstats_pop1_new = {}
    for i,ss_name in enumerate(sumstats_pop1_names):
        ss_df = sumstats_pop1[ss_name].copy()
        ss_df['BETA_tram'] = tram_res_beta[:,i]
        ss_df['SE_tram'] = tram_res_beta_se[:,i]
        ss_df['Z_tram'] = tram_res_beta[:,i]/tram_res_beta_se[:,i]
        ss_df = ss_df.loc[ss_df['CHR'].isin(chrs)]
        ss_df['P_tram'] = calculate_p(ss_df['Z_tram'].values)
        conversion_factor_col = np.reciprocal(np.sqrt(2.0 * ss_df['FRQ'] * (1.0 - ss_df['FRQ'])))
        ss_df['BETA_tram'] = ss_df['BETA_tram']*conversion_factor_col
        ss_df['SE_tram'] = ss_df['SE_tram']*conversion_factor_col
        ss_df = ss_df[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA_tram','SE_tram','Z_tram','P_tram','N']]
        ss_df.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N']
        sumstats_pop1_new[ss_name] = ss_df
    sumstats_pop2_new = {}
    for i,ss_name in enumerate(sumstats_pop2_names):
        ss_df = sumstats_pop2[ss_name].copy()
        ss_df['BETA_tram'] = tram_res_beta[:,i+P1n]
        ss_df['SE_tram'] = tram_res_beta_se[:,i+P1n]
        ss_df['Z_tram'] = tram_res_beta[:,i+P1n]/tram_res_beta_se[:,i+P1n]
        ss_df = ss_df.loc[ss_df['CHR'].isin(chrs)]
        ss_df['P_tram'] = calculate_p(ss_df['Z_tram'].values)
        conversion_factor_col = np.reciprocal(np.sqrt(2.0 * ss_df['FRQ'] * (1.0 - ss_df['FRQ'])))
        ss_df['BETA_tram'] = ss_df['BETA_tram']*conversion_factor_col
        ss_df['SE_tram'] = ss_df['SE_tram']*conversion_factor_col
        ss_df = ss_df[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA_tram','SE_tram','Z_tram','P_tram','N']]
        ss_df.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N']
        sumstats_pop2_new[ss_name] = ss_df

    # Run LD score regressions
    logger.info("Running LD Score regression for the whole genome after Meta-GWAS.")
    Omega_global_new, Sigma_LD_new = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1_new,sumstats_pop2_names,\
        sumstats_pop2_new,df_ldsc_pop1.loc[df_ldsc_pop1['CHR'].isin(chrs),['base']].values,\
            df_ldsc_pop2.loc[df_ldsc_pop2['CHR'].isin(chrs),['base']].values,\
                df_ldsc_te.loc[df_ldsc_te['CHR'].isin(chrs),['base']].values)
    Omega_global_new = Omega_global_new.squeeze()
    logger.info("Regression coefficients (LD):\n%s", Omega_global_new)
    logger.info("Regression coefficients (Intercept):\n%s", Sigma_LD_new)

    logger.info("Estimating the effective sample size and save output to disk")
    for i,ss_name in enumerate(sumstats_pop1_names):
        logger.info("Writing %s (Population 1) to disk", ss_name)
        ss_df_new = sumstats_pop1_new[ss_name]
        ss_df = sumstats_pop1[ss_name]
        Neff_f = ((ss_df_new['Z']**2).mean() - Sigma_LD_new[i,i])/((ss_df['Z']**2).mean() - Sigma_LD[i,i])
        if Neff_f<1:
            Neff_f = ((ss_df_new['Z']**2).mean() - 1)/((ss_df['Z']**2).mean() - Sigma_LD[i,i])
        ss_df_new['N_eff'] = Neff_f*ss_df_new['N'].values
        ss_df_new.to_csv(args.out+'_TRAM_pop1_'+ss_name+'.txt',sep='\t',index=None)

    for i,ss_name in enumerate(sumstats_pop2_names):
        logger.info("Writing %s (Population 2) to disk", ss_name)
        ss_df_new = sumstats_pop2_new[ss_name]
        ss_df = sumstats_pop2[ss_name]
        Neff_f = ((ss_df_new['Z']**2).mean() - Sigma_LD_new[i+P1n,i+P1n])/((ss_df['Z']**2).mean() - Sigma_LD[i+P1n,i+P1n])
        if Neff_f<1:
            Neff_f = ((ss_df_new['Z']**2).mean() - 1)/((ss_df['Z']**2).mean() - Sigma_LD[i+P1n,i+P1n])
        ss_df_new['N_eff'] = Neff_f*ss_df_new['N'].values
        ss_df_new.to_csv(args.out+'_TRAM_pop2_'+ss_name+'.txt',sep='\t',index=None)
    logger.info("\nExecution complete")
