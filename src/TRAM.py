import pandas as pd
import numpy as np
import sys
import argparse
import os
from utils import *

__version__ = '1.0.0'
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
    python <install path>/src/TRAM.py \n
        --out test \n
        --sumstats-popu1 ss_file1,ss_name1 \n
        --sumstats-popu2 ss_file2,ss_name2 \n
        --ldscores ldsc_annot_1mb_TGP_hm3_chr@_std \n
        --annot annot_chr@.gz \
        --out-harmonized \n
        
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
    parser.add_argument('--annot', type=str,
       help='specifies the annotation files used in S-LDXR, \
       If the filename prefix contains the symbol @, TRAM will replace the @ symbol \
            with chromosome number',required=True)
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
            sep='\t',compression='infer',dtype={'CHR':str}).drop_duplicates('SNP').dropna()
        logger.info("Sumstats {} shape: {}x{}".format(ssfile.split(',')[1],
            sumstats_pop1[ssfile.split(',')[1]].shape[0],sumstats_pop1[ssfile.split(',')[1]].shape[1])) 
    sumstats_pop2 = {}
    sumstats_pop2_names = []
    for ssfile in args.sumstats_popu2:
        sumstats_pop2_names.append(ssfile.split(',')[1])
        sumstats_pop2[ssfile.split(',')[1]] = pd.read_csv(ssfile.split(',')[0],
            sep='\t',compression='infer',dtype={'CHR':str}).drop_duplicates('SNP').dropna()
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
        sumstats_pop1[k] = v.loc[v['SNP'].isin(snps_common)].reset_index(drop=True)
    for k,v in sumstats_pop2.items():
        sumstats_pop2[k] = v.loc[v['SNP'].isin(snps_common)].reset_index(drop=True)

    # Removing ambiguous SNPs
    logger.info('Removing ambiguous SNPs...')
    for k,v in sumstats_pop1.items():
        v['A1_int'] = v['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
        v['A2_int'] = v['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
        sumstats_pop1[k] = v
    for k,v in sumstats_pop2.items():
        v['A1_int'] = v['A1'].map(lambda x: 1 if x in ['A','T'] else 2)
        v['A2_int'] = v['A2'].map(lambda x: 1 if x in ['A','T'] else 2)
        sumstats_pop2[k] = v
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
    allele_ref = sumstats_pop1[sumstats_pop1_names[0]][['SNP','A1','A2']]
    allele_ref.columns = ['SNP','A1_ref','A2_ref']
    for k,v in sumstats_pop1.items():
        v = v.merge(allele_ref,left_on='SNP',right_on='SNP',how='inner')
        index_flip = v['A1'].values!=v['A1_ref'].values
        logger.info('{} SNPs need to be flipped for sumstats: {}'.format(index_flip.sum(),k))
        v.loc[index_flip,'BETA'] = -v.loc[index_flip,'BETA']
        v.loc[index_flip,'Z'] = -v.loc[index_flip,'Z']
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

    # Run LD score regressions for annotation
    logger.info("Running LD Score regression for local regions/functional annotations chromosome by chromosome.")
    if args.allowed_chr_values is None:
        chrs = [str(_) for _ in range(1,23)]
    else:
        chrs = args.allowed_chr_values

    # beta, se, N_eff
    tram_res_beta = np.zeros((df_ldsc_pop1.shape[0],P1n+P2n))
    tram_res_beta_se = np.zeros((df_ldsc_pop1.shape[0],P1n+P2n))
    tram_res_Neff = np.zeros((df_ldsc_pop1.shape[0],P1n+P2n))
    res_omegas_local = {}

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
            se_array[:,i+P2n] = v['SE'].values

        name_windows = list(ldscore1_c.columns.values[4:])
        n_window = len(name_windows)
        ldscore1_ck = ldscore1_c[['base']+name_windows].values
        ldscore2_ck = ldscore2_c[['base']+name_windows].values
        ldscorete_ck = ldscorete_c[['base']+name_windows].values
        if args.local_anno == 'yes':
            logger.info("find {} local regions in chromosome {}".format(len(name_windows),c))
            omegas_local = np.zeros((1+len(name_windows),P1n+P2n,P1n+P2n))
            omegas_local_baseLD = np.zeros((len(name_windows),P1n+P2n,P1n+P2n))
            for i,name_window in enumerate(name_windows):
                omegas_local_k, _ = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1_c,sumstats_pop2_names,\
                    sumstats_pop2_c,ldscore1_ck[:,[0,i+1]],ldscore2_ck[:,[0,i+1]],ldscorete_ck[:,[0,i+1]],Sigma_LD)
                omegas_local[i+1] = make_positive_definite(omegas_local_k.sum(axis=0))
                omegas_local_baseLD[i] = omegas_local_k[0]
            omegas_local_baseLD_mean = make_positive_definite(omegas_local_baseLD.mean(axis=0))
            omegas_local[0] = 2*omegas_local_baseLD_mean
            omegas_local -= omegas_local_baseLD_mean
        else:
            # 1+K x P x P
            logger.info("find {} functional annotations in chromosome {}".format(len(name_windows),c))
            omegas_local, _ = run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1_c,sumstats_pop2_names,\
                sumstats_pop2_c,ldscore1_ck,ldscore2_ck,ldscorete_ck,Sigma_LD)
        logger.info("Local regions/functional annotation regression coefficients (LD):\n%s \n...", omegas_local[:2])
        res_omegas_local[c] = omegas_local

        logger.info("Creating omega and sigma matrices for each SNP.")
        omega_ldscore_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))
        sigma_snps = np.zeros((snp_num_c,P1n+P2n,P1n+P2n))
        for j in range(snp_num_c):
            sigma_snps[j] = np.diag(se_array[j])@Sigma_LD@np.diag(se_array[j])
            ld_matrx_j = np.zeros_like(omegas_local)
            for k in range(n_window+1):
                ld_matrx_j[k][:P1n,:P1n] = ldscore1_ck[j,k]
                ld_matrx_j[k][:P1n,P1n:] = ld_matrx_j[k][P1n:,:P1n] = ldscorete_ck[j,k]
                ld_matrx_j[k][P1n:,P1n:] = ldscore2_ck[j,k]
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
        logger.info("Run the TRAM method\n")
        new_betas, new_beta_ses = run_tram_method(beta_array,omega_ldscore_snps,sigma_snps)
        new_Neff = np.zeros((snp_num_c,P1n+P2n))
        for i,ss_name in enumerate(sumstats_pop1_names):
            new_Neff[:,i] = calculate_n_eff(i,sumstats_pop1_c[ss_name]['N'].values,sigma_snps,new_beta_ses[:,i])
        for i,ss_name in enumerate(sumstats_pop2_names):
            new_Neff[:,P1n+i] = calculate_n_eff(P1n+i,sumstats_pop2_c[ss_name]['N'].values,sigma_snps,new_beta_ses[:,P1n+i])
        tram_res_beta[chr_snp_idx] = new_betas
        tram_res_beta_se[chr_snp_idx] = new_beta_ses
        tram_res_Neff[chr_snp_idx] = new_Neff

    if args.out_reg_coef == 'yes':
        logger.info("Saving local regions/functional annotation regression coefficients to file %s",args.out+'_TRAM_reg_coefs.npy')
        np.save(args.out+'_TRAM_reg_coefs.npy',res_omegas_local)

    logger.info("Preparing results for output")
    for i,ss_name in enumerate(sumstats_pop1_names):
        logger.info("Writing %s (Population 1) to disk", ss_name)
        ss_df = sumstats_pop1[ss_name]
        ss_df['BETA_tram'] = tram_res_beta[:,i]
        ss_df['SE_tram'] = tram_res_beta_se[:,i]
        ss_df['Z_tram'] = tram_res_beta[:,i]/tram_res_beta_se[:,i]
        ss_df['P_tram'] = calculate_p(ss_df['Z_tram'].values)
        ss_df['N_eff'] = tram_res_Neff[:,i]
        conversion_factor_col = np.reciprocal(np.sqrt(2.0 * ss_df['FRQ'] * (1.0 - ss_df['FRQ'])))
        ss_df['BETA_tram'] = ss_df['BETA_tram']*conversion_factor_col
        ss_df['SE_tram'] = ss_df['SE_tram']*conversion_factor_col
        ss_df = ss_df[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA_tram','SE_tram','Z_tram','P_tram','N_eff','N']]
        ss_df.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N_eff','N_origin']
        ss_df.to_csv(args.out+'_TRAM_pop1_'+ss_name+'.txt',sep='\t',index=None)
    for i,ss_name in enumerate(sumstats_pop2_names):
        logger.info("Writing %s (Population 2) to disk", ss_name)
        ss_df = sumstats_pop2[ss_name]
        ss_df['BETA_tram'] = tram_res_beta[:,i+P1n]
        ss_df['SE_tram'] = tram_res_beta_se[:,i+P1n]
        ss_df['Z_tram'] = tram_res_beta[:,i+P1n]/tram_res_beta_se[:,i+P1n]
        ss_df['P_tram'] = calculate_p(ss_df['Z_tram'].values)
        ss_df['N_eff'] = tram_res_Neff[:,i+P1n]
        conversion_factor_col = np.reciprocal(np.sqrt(2.0 * ss_df['FRQ'] * (1.0 - ss_df['FRQ'])))
        ss_df['BETA_tram'] = ss_df['BETA_tram']*conversion_factor_col
        ss_df['SE_tram'] = ss_df['SE_tram']*conversion_factor_col
        ss_df = ss_df[['CHR','BP','SNP','A1_ref','A2_ref','FRQ','BETA_tram','SE_tram','Z_tram','P_tram','N_eff','N']]
        ss_df.columns = ['CHR','BP','SNP','A1','A2','FRQ','BETA','SE','Z','P','N_eff','N_origin']
        ss_df.to_csv(args.out+'_TRAM_pop2_'+ss_name+'.txt',sep='\t',index=None)
    logging.info("\nExecution complete")