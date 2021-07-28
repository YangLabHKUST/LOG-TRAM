import pandas as pd
import numpy as np
import logging
import sys
import copy
from pandas_plink import read_plink1_bin
import os
from scipy import linalg
import scipy.stats as st
from scipy.stats import norm
import gc
import warnings
warnings.filterwarnings('ignore') 


##
# Data loading and preprocessing
##
def configure_logging(logger_name): 
    LOG_LEVEL = logging.INFO
    log_filename = logger_name+'.log'
    importer_logger = logging.getLogger('importer_logger')
    importer_logger.setLevel(LOG_LEVEL)
    formatter = logging.Formatter('%(asctime)s : %(levelname)s : %(message)s')

    fh = logging.FileHandler(filename=log_filename)
    fh.setLevel(LOG_LEVEL)
    fh.setFormatter(formatter)
    importer_logger.addHandler(fh)

    sh = logging.StreamHandler(sys.stdout)
    sh.setLevel(LOG_LEVEL)
    sh.setFormatter(formatter)
    importer_logger.addHandler(sh)
    return importer_logger


def ReadPlink(plink_file,bim,dtype=np.float32):
    Genotype = read_plink1_bin(plink_file+".bed", plink_file+".bim", plink_file+".fam", verbose=False)
    Genotype = Genotype.where(Genotype.snp.isin(Genotype.snp.values[bim['index'].values]), drop=True)
    Genotype = Genotype.astype(np.int8)
    G_geno = Genotype.values
    G_geno[np.isnan(G_geno)] = 2
    G_geno = 2-G_geno
    return G_geno.astype(dtype)


def ScaleGenotype(G_geno,dtype=np.float32):
    # scale to mean 0, variance 1/p
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    G_geno = ((G_geno-G_geno_mean)/G_geno_sd)/np.sqrt(G_geno.shape[1])
    return G_geno.astype(dtype), G_geno_mean/2, G_geno_sd


def GenotypeMoment(G_geno):
    G_geno_mean = G_geno.mean(axis=0)
    G_geno_sd = G_geno.var(axis=0)**.5
    return G_geno_mean, G_geno_sd


def PLINK2MAMA(file,frq):
    # file: plink2 gwas output
    # fre: plink --freq output (for extract maf)
    df = pd.read_csv(file,sep='\t',compression='infer')
    frq = pd.read_csv(frq,delim_whitespace=True)
    df = df.merge(frq[['SNP','MAF']],left_on='ID',right_on='SNP')
    df['z'] = df.apply(lambda x:x['BETA']/x['SE'], axis=1)
    df = df[['ID','#CHROM','POS','A1','REF','MAF','BETA','SE','z','P','OBS_CT']]
    df.columns = ['SNP','CHR','BP','A1','A2','FRQ','BETA','SE','Z','P','N']
    return df


def BOLT2MAMA(file,n):
    # file: bolt gwas output
    # n: GWAS sample size
    df = pd.read_csv(file,sep='\t',compression='infer')
    df['n'] = n-df['F_MISS']*n
    df['n'] = df['n'].astype(int)
    df['z'] = df.apply(lambda x:x['BETA']/x['SE'], axis=1)
    try:
        df = df[['SNP','CHR','BP','ALLELE1','ALLELE0','A1FREQ','BETA','SE','z','P_BOLT_LMM','n']]
    except:
        df = df[['SNP','CHR','BP','ALLELE1','ALLELE0','A1FREQ','BETA','SE','z','P_BOLT_LMM_INF','n']]
    df.columns = ['SNP','CHR','BP','A1','A2','FRQ','BETA','SE','Z','P','N']
    return df


##
# refer to LDXR
##
def get_coef_raw(ldscore_mat, zsc1, zsc2, max_int):
    """
    Obtain coefficient for heritability / genetic covariance
    """
   
    # get dimension
    nsnp,ncoef = ldscore_mat.shape[0],ldscore_mat.shape[1]

    # set up variables
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    score_sum = np.sum(ldscore_mat[:,:-1], axis=1)
    
    # get averages
    mean_ldscore = np.mean(score_sum)
    mean_zprod = np.mean(zprod)
    mean_nprod = np.mean(nprod)

    # estimate intercept
    zprod_sorted = np.sort(zprod)
    idx = int(float(nsnp)*0.95)
    intercept = np.mean(zprod_sorted[0:idx])
    if intercept > max_int:
        intercept = max_int

    # get raw coef
    coef = (mean_zprod-intercept) / (mean_nprod*mean_ldscore)

    return coef, intercept


def get_pred(coef, x, n1, n2, intercept):
    
    nprod = np.sqrt(n1*n2)
    score_sum = np.sum(x[:,:-1], axis=1)
    pred = coef*score_sum*nprod + intercept

    return pred.astype(np.float32)


def update_weights(w1, w2, wx, pred1, pred2, predx):
   
    var1 = 2.0*np.square(pred1)
    var2 = 2.0*np.square(pred2)
    var12 = np.sqrt(var1*var2)/2.0 + np.square(predx)

    new_w1 = w1 / var1
    new_w2 = w2 / var2
    new_wx = wx / var12

    return new_w1, new_w2, new_wx


def regression(x, y, constrain_intercept,subtract):
    """
    Perform least square regression with block jack knife
    """

    # perform regression
    xtx = np.dot(x.T, x)
    xty = np.dot(x.T, y)
    ncoef = xtx.shape[0]

    # obtain coefficient
    coef = np.zeros(ncoef, dtype=np.float32)
    if constrain_intercept:
        coef[:-1] = np.linalg.solve(xtx[:-1,:-1], xty[:-1])
        coef[-1] = subtract
    else:
        coef = np.linalg.solve(xtx, xty)

    return coef


def get_coef(score, zsc1, zsc2, w, constrain_intercept, subtract):
    """
    Obtain coefficient for heritability / genetic covariance
    """
    
    # set up the regression
    nsnp = score.shape[0]
    nprod = np.sqrt((zsc1['N'].values)*(zsc2['N'].values))
    zprod = (zsc1['Z'].values)*(zsc2['Z'].values)
    if constrain_intercept:
        zprod -= subtract
    score[:,:-1] = score[:,:-1] * nprod[:,np.newaxis]
    
    # scale the matrix to improve matrix condition
    nbar = np.mean(nprod)
    score[:,:-1] /= nbar
    
    # apply the weight
    score_w = score * np.sqrt(w[:,np.newaxis])
    zprod_w = zprod * np.sqrt(w)
    
    # run regression
    coef = regression(score_w, zprod_w, constrain_intercept, subtract)
    coef[:-1] /= nbar
    
    return coef
  

def estimate_h2(sumstat1,ldscore1,reg_w1=1,constrain_intercept=False,int1=1):
    reg_w1[reg_w1<0] = 0.01
    # tau: [coeficient, intercept]
    ldscore1 = np.hstack((ldscore1,np.ones((ldscore1.shape[0],1))))
    tau1,int1_ = get_coef_raw(ldscore1,sumstat1,sumstat1,1)
    n1 = sumstat1['N'].values
    pred1 = get_pred(tau1, ldscore1, n1, n1, int1_)
    var1 = 2.0*np.square(pred1)
    weight1_ = reg_w1 / var1
    tau1 = get_coef(ldscore1, sumstat1, sumstat1,
       weight1_, constrain_intercept, int1)
    return tau1


def estimate_gc(sumstat1,sumstat2,ldscore1,ldscore2,ldscorex,reg_w1=1,reg_w2=1,constrain_intercept=False,int1=1,int2=1,intx=0):
    reg_w1[reg_w1<0] = 0.01
    reg_w2[reg_w2<0] = 0.01
    reg_wx = np.sqrt(reg_w1*reg_w2)

    ldscore1 = np.hstack((ldscore1,np.ones((ldscore1.shape[0],1))))
    ldscore2 = np.hstack((ldscore2,np.ones((ldscore2.shape[0],1))))
    ldscorex = np.hstack((ldscorex,np.ones((ldscorex.shape[0],1))))

    # get initial estimate of coefs for constructing weights
    tau1,int1_ = get_coef_raw(ldscore1,sumstat1,sumstat1,1)
    tau2,int2_ = get_coef_raw(ldscore2,sumstat2,sumstat2,1)
    theta,intx_ = get_coef_raw(ldscorex,sumstat1,sumstat2,0)

    # predict chi-square and prodcut 
    n1,n2 = sumstat1['N'].values,sumstat2['N'].values
    pred1 = get_pred(tau1, ldscore1, n1, n1, int1_)
    pred2 = get_pred(tau2, ldscore2, n2, n2, int2_)
    predx = get_pred(theta, ldscorex, n1, n2, intx_)
    
    # update weight
    weight1_, weight2_, weightx_ = update_weights(reg_w1, reg_w2, reg_wx, pred1, pred2, predx)

    # get regression coefficients
    tau1 = get_coef(ldscore1, sumstat1, sumstat1,
       weight1_, constrain_intercept, int1)
    tau2 = get_coef(ldscore2, sumstat2, sumstat2,
        weight2_, constrain_intercept, int2)
    theta = get_coef(ldscorex, sumstat1, sumstat2,
        weightx_, constrain_intercept, intx)

    return tau1,tau2,theta


def ldscore_regression_h2(df_eas,ldscore,intercept=None):
    if intercept is None:
        idx1 = (df_eas['Z']**2<30).values
        fit_step1 = estimate_h2(df_eas.loc[idx1],ldscore[idx1],
                reg_w1 = 1/ldscore[idx1,0],constrain_intercept=False)
        fit_step2 = estimate_h2(df_eas,ldscore,
                reg_w1 = 1/ldscore[:,0],constrain_intercept=True,int1=fit_step1[-1])
    else:
        fit_step2 = estimate_h2(df_eas,ldscore,
                reg_w1 = 1/ldscore[:,0],constrain_intercept=True,int1=intercept)
    fit_step2[0] = fit_step2[0] if fit_step2[0]>0 else 1e-12 

    return fit_step2


def ldscore_regression_gc(df_eas,df_eur,ldscore1,ldscore2,ldscorete,intercept1=None,intercept2=None,interceptx=None):
    idx1 = (df_eas['Z']**2<30).values & (df_eur['Z']**2<30).values
    fit_step1 = estimate_gc(df_eas.loc[idx1],df_eur.loc[idx1],ldscore1[idx1],ldscore2[idx1],ldscorete[idx1],
            reg_w1 = 1/ldscore1[idx1,0],reg_w2 = 1/ldscore2[idx1,0],constrain_intercept=False)
    if intercept1 is None:
        intercept1 = fit_step1[0][-1]
    if intercept2 is None:
        intercept2 = fit_step1[1][-1]
    if interceptx is None:
        interceptx = fit_step1[2][-1]
    fit_step2 = estimate_gc(df_eas,df_eur,ldscore1,ldscore2,ldscorete,
            reg_w1 = 1/ldscore1[:,0],reg_w2 = 1/ldscore2[:,0],constrain_intercept=True,
            int1=intercept1,int2=intercept2,intx=interceptx)
    
    h1_snp = fit_step2[0][0]
    h2_snp = fit_step2[1][0]
    
    h1_snp = h1_snp if h1_snp>0 else 1e-12 
    h2_snp = h2_snp if h2_snp>0 else 1e-12 
    gc_avg = fit_step2[2][0]/np.sqrt(h1_snp*h2_snp)
    gc_avg = gc_avg if gc_avg<1 else 0.95
    gc_avg = gc_avg if gc_avg>-1 else -0.95   
    h12_snp = np.sqrt(h1_snp*h2_snp)*gc_avg
    fit_step2[2][0] = h12_snp
    return fit_step2[2]


def run_ldscore_regressions(sumstats_pop1_names,sumstats_pop1,sumstats_pop2_names,sumstats_pop2,ldscore1,ldscore2,ldscorete,intercept=None):
    '''
    :return: A tuple holding regression coefficient matrices (ldscores, intercept)
    '''
    P1, P2 = len(sumstats_pop1_names), len(sumstats_pop2_names)
    result_coefs = np.zeros((ldscore1.shape[1]+1, P1+P2, P1+P2))
    if intercept is not None:
        if intercept.shape[0] != intercept.shape[1] or intercept.shape[0] != P1+P2:
            raise ValueError('user definded intercept must be a sqaure matrix with dimention = number of susmstats')

    # Calculate diagonal coefficients first
    for p,ss_name in enumerate(sumstats_pop1_names):
        if intercept is None:
            result_coefs[:,p,p] = ldscore_regression_h2(sumstats_pop1[ss_name],ldscore1)
        else:
            result_coefs[:,p,p] = ldscore_regression_h2(sumstats_pop1[ss_name],ldscore1,intercept[p,p])
    for p,ss_name in enumerate(sumstats_pop2_names):
        if intercept is None:
            result_coefs[:,p+P1,p+P1] = ldscore_regression_h2(sumstats_pop2[ss_name],ldscore2)
        else:
            result_coefs[:,p+P1,p+P1] = ldscore_regression_h2(sumstats_pop2[ss_name],ldscore2,intercept[p+P1,p+P1])
    

    # Calculate each off-diagonal element
    for i in range(P1):
        for j in range(i+1,P1):
            if intercept is None:
                result_coefs[:,i,j] = ldscore_regression_gc(sumstats_pop1[sumstats_pop1_names[i]],\
                    sumstats_pop1[sumstats_pop1_names[j]],ldscore1,ldscore1,ldscore1)
            else:
                result_coefs[:,i,j] = ldscore_regression_gc(sumstats_pop1[sumstats_pop1_names[i]],\
                    sumstats_pop1[sumstats_pop1_names[j]],ldscore1,ldscore1,ldscore1,\
                        intercept[i,i],intercept[j,j],intercept[i,j])
            result_coefs[:,j,i] = result_coefs[:,i,j] 
    for i in range(P1):
        for j in range(P2):
            if intercept is None:
                result_coefs[:,i,P1+j] = ldscore_regression_gc(sumstats_pop1[sumstats_pop1_names[i]],\
                    sumstats_pop2[sumstats_pop2_names[j]],ldscore1,ldscore2,ldscorete)
            else:
                result_coefs[:,i,P1+j] = ldscore_regression_gc(sumstats_pop1[sumstats_pop1_names[i]],\
                    sumstats_pop2[sumstats_pop2_names[j]],ldscore1,ldscore2,ldscorete,\
                        intercept[i,i],intercept[P1+j,P1+j],intercept[i,P1+j])
            result_coefs[:,P1+j,i] = result_coefs[:,i,P1+j]
    for i in range(P2):
        for j in range(i+1,P2):
            if intercept is None:
                result_coefs[:,P1+i,P1+j] = ldscore_regression_gc(sumstats_pop2[sumstats_pop2_names[i]],\
                    sumstats_pop2[sumstats_pop2_names[j]],ldscore2,ldscore2,ldscore2)
            else:
                result_coefs[:,P1+i,P1+j] = ldscore_regression_gc(sumstats_pop2[sumstats_pop2_names[i]],\
                    sumstats_pop2[sumstats_pop2_names[j]],ldscore2,ldscore2,ldscore2,\
                        intercept[P1+i,P1+i],intercept[P1+j,P1+j],intercept[P1+i,P1+j])
            result_coefs[:,P1+j,P1+i] = result_coefs[:,P1+i,P1+j]
    return result_coefs[:-1], result_coefs[-1]    


def make_positive_definite(A):
    d = A.shape[0]
    for i in range(d):
        A[i,i] = A[i,i] if A[i,i]>0 else 1e-12 
    for i in range(d):
        for j in range(i+1,d):
            if A[i,i]<=1e-12 or A[j,j]<=1e-12:
                A[i,j] = A[j,i] = 0
                continue
            if A[i,j] > 0.95*(A[i,i]*A[j,j])**.5:
                A[i,j] = A[j,i] = 0.95*(A[i,i]*A[j,j])**.5
            if A[i,j] < -0.95*(A[i,i]*A[j,j])**.5:
                A[i,j] = A[j,i] = -0.95*(A[i,i]*A[j,j])**.5
    return A


def get_enrichment(coef, annot_mat):
    
    """
    Estimate enrichment of local regions/functional annotation
    """
       
    # exclude coef for the intercept 
    annot_cov = np.dot(annot_mat.T, annot_mat)  
    annot_est = np.dot(annot_cov, coef)

    annot_cov = annot_cov.astype(np.float64)
    annot_est = annot_est.astype(np.float64)
    
    annot_nsnp = np.sum(annot_mat, axis=0)
    # get dimension
    tot_nsnp = np.float64(annot_nsnp[0])
    annot_en = annot_est*tot_nsnp / (annot_est[0] * annot_nsnp)

    return annot_en[1:]


##
# Association mapping
##
def calculate_p(z_scores: np.array) -> np.array:
    # Calculate constants used in determination of P values
    ln = np.log  # pylint: disable=invalid-name
    LN_2 = ln(2.0)
    RECIP_LN_10 = np.reciprocal(ln(10.0))
    """
    Function that calculates P for the TRAM results
    :param z_scores: Z scores
    :return: P values 
             (as strings, to allow for very large negative exponents)
    """
    # Since P = 2 * normal_cdf(-|Z|), P = e ^ (log_normal_cdf(-|Z|) + ln 2)
    # This can be changed to base 10 as P = 10 ^ ((log_normal_cdf(-|Z|) + ln 2) / ln 10)
    log_10_p = RECIP_LN_10 * (norm.logcdf(-np.abs(z_scores)) + LN_2)

    # Break up the log based 10 of P values into the integer and fractional part
    # To handle the case of Z = 0 (and not result in "10e-1"), set initial values to (-1.0, 1.0)
    frac_part, int_part = np.full_like(z_scores, -1.0), np.full_like(z_scores, 1.0)
    np.modf(log_10_p, out=(frac_part, int_part), where=(z_scores != 0.0))

    # Construct strings for the P values
    # 1) Add one to the fractional part to ensure that the result mantissa is between 1 and 10
    # 2) Subtract one from the integer part to compensate and keep the overall value correct
    result = np.char.add(np.char.add(np.power(10.0, (frac_part + 1.0)).astype(str), 'e'),
                         (int_part - 1).astype(int).astype(str))

    return result


def MTAG(df_res,beta_array,se_array,ldscore1,ldscore2,ldscorete,omega,sigma_ld,verbose=False):
    p = beta_array.shape[0]
    blue_estimate_t = np.zeros(p)
    blue_estimate_s = np.zeros(p)
    blue_estimate_t_sd = np.zeros(p)
    blue_estimate_s_sd = np.zeros(p)
    pval_t = np.zeros(p)
    pval_s = np.zeros(p)
    #beta_combine = np.vstack((df_eas['BETA'],df_eur['BETA'])).T
    #se_combine = np.vstack((df_eas['SE'],df_eur['SE'])).T
    #omega_mtag = np.zeros((beta_combine.shape[0],2,2))
    #for i in range(beta_combine.shape[0]):
    #    W_i = np.diag(se_combine[i])
    #    omega_mtag[i] = beta_combine[i].reshape(2,1)@beta_combine[i].reshape(1,2)-W_i@sigma_ld@W_i
    #omega_j = omega_mtag.mean(axis=0)
    #if not np.all(linalg.eigvalsh(omega_j) >= 0.0):
    omega_j = omega*np.array([[np.median(ldscore1),np.median(ldscorete)],[np.median(ldscorete),np.median(ldscore2)]])
    if not np.all(linalg.eigvalsh(omega_j) >= 0.0):
        raise Exception('Omega is not semi-positive')
    for j in range(p):
        if j%10000==0 and j!=0 and verbose:
            print(j,end=',')
        beta_j_hat = beta_array[j]
        # construct Omega and Sigma
        W_j = np.diag(se_array[j])
        sigma_j = W_j@sigma_ld@W_j
        # Check omega and sigma for validity based on positive (semi-)definiteness
        try:
            linalg.cholesky(sigma_j)
        except linalg.LinAlgError:
            blue_estimate_t[j] = np.nan
            continue
        omega_j_t_vector_norm = omega_j[:,0].reshape(-1,1)/omega_j[0,0]
        omega_j_s_vector_norm = omega_j[:,1].reshape(-1,1)/omega_j[1,1]
        lambda_t = linalg.inv(omega_j+sigma_j-omega_j[:,0].reshape(-1,1)@omega_j[:,0].reshape(-1,1).T/omega_j[0,0])
        lambda_s = linalg.inv(omega_j+sigma_j-omega_j[:,1].reshape(-1,1)@omega_j[:,1].reshape(-1,1).T/omega_j[1,1])
        blue_estimate_t[j] = omega_j_t_vector_norm.T@lambda_t@beta_j_hat/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)
        blue_estimate_s[j] = omega_j_s_vector_norm.T@lambda_s@beta_j_hat/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)
        blue_estimate_t_sd[j] = 1/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)[0,0]**.5
        blue_estimate_s_sd[j] = 1/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)[0,0]**.5
        pval_t[j] = st.norm.sf(abs(blue_estimate_t[j]/blue_estimate_t_sd[j]))*2
        pval_s[j] = st.norm.sf(abs(blue_estimate_s[j]/blue_estimate_s_sd[j]))*2
        
    #df_res = copy.deepcopy(df_eas[['SNP','CHR','BP','A1','A2']])       
    df_res['beta_EAS_homo'] = blue_estimate_t
    df_res['se_EAS_homo'] = blue_estimate_t_sd
    df_res['pval_EAS_homo'] = pval_t
    df_res['beta_EUR_homo'] = blue_estimate_s
    df_res['se_EUR_homo'] = blue_estimate_s_sd
    df_res['pval_EUR_homo'] = pval_s   
    df_res = df_res.dropna()
    
    return df_res


def MAMA(df_res,beta_array,se_array,ldscore1,ldscore2,ldscorete,omega,sigma_ld,verbose=False):
    p = df_res.shape[0]
    blue_estimate_t = np.zeros(p)
    blue_estimate_s = np.zeros(p)
    blue_estimate_t_sd = np.zeros(p)
    blue_estimate_s_sd = np.zeros(p)
    pval_t = np.zeros(p)
    pval_s = np.zeros(p)
    for j in range(p):
        if j%10000==0 and j!=0 and verbose:
            print(j,end=',')
        beta_j_hat = beta_array[j]
        # construct Omega and Sigma
        omega_j = omega*np.array([[ldscore1[j],ldscorete[j]],[ldscorete[j],ldscore2[j]]])
        W_j = np.diag(se_array[j])
        sigma_j = W_j@sigma_ld@W_j
        # Check omega and sigma for validity based on positive (semi-)definiteness
        if not np.all(linalg.eigvalsh(omega_j) >= 0.0):
            blue_estimate_t[j] = np.nan
            continue
        try:
            linalg.cholesky(sigma_j)
        except linalg.LinAlgError:
            blue_estimate_t[j] = np.nan
            continue
        omega_j_t_vector_norm = omega_j[:,0].reshape(-1,1)/omega_j[0,0]
        omega_j_s_vector_norm = omega_j[:,1].reshape(-1,1)/omega_j[1,1]
        lambda_t = linalg.inv(omega_j+sigma_j-omega_j[:,0].reshape(-1,1)@omega_j[:,0].reshape(-1,1).T/omega_j[0,0])
        lambda_s = linalg.inv(omega_j+sigma_j-omega_j[:,1].reshape(-1,1)@omega_j[:,1].reshape(-1,1).T/omega_j[1,1])
        blue_estimate_t[j] = omega_j_t_vector_norm.T@lambda_t@beta_j_hat/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)
        blue_estimate_s[j] = omega_j_s_vector_norm.T@lambda_s@beta_j_hat/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)
        blue_estimate_t_sd[j] = 1/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)[0,0]**.5
        blue_estimate_s_sd[j] = 1/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)[0,0]**.5
        pval_t[j] = st.norm.sf(abs(blue_estimate_t[j]/blue_estimate_t_sd[j]))*2
        pval_s[j] = st.norm.sf(abs(blue_estimate_s[j]/blue_estimate_s_sd[j]))*2
        
    #df_res = copy.deepcopy(df_eas[['SNP','CHR','BP','A1','A2']])       
    df_res['beta_EAS_homo'] = blue_estimate_t
    df_res['se_EAS_homo'] = blue_estimate_t_sd
    df_res['pval_EAS_homo'] = pval_t
    df_res['beta_EUR_homo'] = blue_estimate_s
    df_res['se_EUR_homo'] = blue_estimate_s_sd
    df_res['pval_EUR_homo'] = pval_s   
    df_res = df_res.dropna()

    return df_res


def AnnoWindowParam(df_res,df_param_c,c):
    '''
    annotate local omega/simga for each SNP in chromosome 'c'
    '''
    for row in df_param_c.iterrows():
        start = int(row[1]['start'])
        end = int(row[1]['end'])
        if pd.isnull(row[1]['h1_perSNP_window']) or pd.isnull(row[1]['h2_perSNP_window']):
            df_res.loc[(df_res['CHR']==c)&(df_res['BP']>=start)&(df_res['BP']<end),
                    ['h1_perSNP_window','h2_perSNP_window','h12_perSNP_window','enrichment_h1','enrichment_h2','enrichment_h12']]= \
            0,0,0,0,0,0
        else:
            h1_window = row[1]['h1_perSNP_window']
            h2_window = row[1]['h2_perSNP_window']
            h12_window = row[1]['h12_perSNP_window']
            df_res.loc[(df_res['CHR']==c)&(df_res['BP']>=start)&(df_res['BP']<end),
                        ['h1_perSNP_window','h2_perSNP_window','h12_perSNP_window','enrichment_h1','enrichment_h2','enrichment_h12']]= \
                        h1_window,h2_window,h12_window,row[1]['enrichment_h1'],row[1]['enrichment_h2'],row[1]['enrichment_h12']
    return df_res  


def TramNaive(df_anno_params,beta_array,se_array,ldscore1,ldscore2,ldscorete,sigma_ld,verbose=False):
    p = df_anno_params.shape[0]
    blue_estimate_t = np.zeros(p)
    blue_estimate_s = np.zeros(p)
    blue_estimate_t_sd = np.zeros(p)
    blue_estimate_s_sd = np.zeros(p)
    pval_t = np.zeros(p)
    pval_s = np.zeros(p)
    h1_windows = df_anno_params['h1_perSNP_window'].values
    h2_windows = df_anno_params['h2_perSNP_window'].values
    h12_windows = df_anno_params['h12_perSNP_window'].values
    for j in range(p):
        if j%10000==0 and j!=0 and verbose:
            print(j,end=',')
        beta_j_hat = beta_array[j]
        # construct Omega and Sigma
        h1_window_j = h1_windows[j]
        h1_window_j = h1_window_j if h1_window_j>0 else 1e-12
        h2_window_j = h2_windows[j]
        h2_window_j = h2_window_j if h2_window_j>0 else 1e-12
        omega_j = np.array([[h1_window_j,h12_windows[j]],
                          [h12_windows[j],h2_window_j]])
        omega_j = omega_j*np.array([[ldscore1[j],ldscorete[j]],[ldscorete[j],ldscore2[j]]])
        W_j = np.diag(se_array[j])
        sigma_j = W_j@sigma_ld@W_j
        # Check omega and sigma for validity based on positive (semi-)definiteness
        if not np.all(linalg.eigvalsh(omega_j) >= 0.0):
            blue_estimate_t[j] = np.nan
            continue
        try:
            linalg.cholesky(sigma_j)
        except linalg.LinAlgError:
            blue_estimate_t[j] = np.nan
            continue
        omega_j_t_vector_norm = omega_j[:,0].reshape(-1,1)/omega_j[0,0]
        omega_j_s_vector_norm = omega_j[:,1].reshape(-1,1)/omega_j[1,1]
        lambda_t = linalg.inv(omega_j+sigma_j-omega_j[:,0].reshape(-1,1)@omega_j[:,0].reshape(-1,1).T/omega_j[0,0])
        lambda_s = linalg.inv(omega_j+sigma_j-omega_j[:,1].reshape(-1,1)@omega_j[:,1].reshape(-1,1).T/omega_j[1,1])
        blue_estimate_t[j] = omega_j_t_vector_norm.T@lambda_t@beta_j_hat/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)
        blue_estimate_s[j] = omega_j_s_vector_norm.T@lambda_s@beta_j_hat/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)
        blue_estimate_t_sd[j] = 1/(omega_j_t_vector_norm.T@lambda_t@omega_j_t_vector_norm)[0,0]**.5
        blue_estimate_s_sd[j] = 1/(omega_j_s_vector_norm.T@lambda_s@omega_j_s_vector_norm)[0,0]**.5
        pval_t[j] = st.norm.sf(abs(blue_estimate_t[j]/blue_estimate_t_sd[j]))*2
        pval_s[j] = st.norm.sf(abs(blue_estimate_s[j]/blue_estimate_s_sd[j]))*2
        
    df_res = copy.deepcopy(df_anno_params[['SNP','CHR','BP','A1','A2']])       
    df_res['beta_EAS_het'] = blue_estimate_t
    df_res['se_EAS_het'] = blue_estimate_t_sd
    df_res['pval_EAS_het'] = pval_t
    df_res['beta_EUR_het'] = blue_estimate_s
    df_res['se_EUR_het'] = blue_estimate_s_sd
    df_res['pval_EUR_het'] = pval_s   
    df_res = df_res.dropna()

    return df_res


def run_tram_method(betas, omega, sigma):
    """
    Runs the core tram method to combine results and generate final, combined summary statistics
    :param betas: MxP matrix of beta values (M = # of SNPs, P = # of ancestries)
    :param omega: MxPxP matrix of omega values (M = # of SNPs, P = # of ancestries)
    :param sigma: MxPxP matrix of sigma values (M = # of SNPs, P = # of ancestries)
    :return: Tuple containing:
                 1) Result ndarray of betas (MxP) where M = SNPs and P = populations
                 2) Result ndarray of beta standard errors (MxP) where M = SNPs and P = populations
    """

    # Get values for M and P (used to keep track of slices / indices / broadcasting)
    M, P, *extra_dimensions = omega.shape  # pylint: disable=unused-variable

    # Create a 3D matrix, M rows of Px1 column vectors with shape (M, P, 1)
    d_indices = np.arange(P)
    omega_diag = omega[:, d_indices, d_indices][:, :, np.newaxis]
    omega_pp_scaled = np.divide(omega, omega_diag)  # Slice rows are Omega'_pjj / omega_pp,j

    # Produce center matrix in steps (product of omega terms, add omega and sigma, then invert)
    center_matrix_inv = -omega_pp_scaled[:, :, :, np.newaxis] * omega[:, :, np.newaxis, :]
    center_matrix_inv += omega[:, np.newaxis, :, :] + sigma[:, np.newaxis, :, :] # Broadcast add
    center_matrix = np.linalg.inv(center_matrix_inv) # Inverts each slice separately
    del center_matrix_inv  # Clean up the inverse matrix to free space
    gc.collect()

    # Calculate (Omega'_p,j/omega_pp,j) * center_matrix
    left_product = np.matmul(omega_pp_scaled[:, :, np.newaxis, :], center_matrix)
    del center_matrix  # Clean up the center matrix to free space
    gc.collect()

    # Calculate denominator (M x P x 1 x 1)
    denom = np.matmul(left_product,
                      np.transpose(omega_pp_scaled[:, :, np.newaxis, :], (0, 1, 3, 2)))
    denom_recip = np.reciprocal(denom)
    denom_recip_view = denom_recip.view()
    denom_recip_view.shape = (M, P)

    # Calculate numerator (M x P x 1 x 1))
    left_product_view = left_product.view().reshape(M, P, P)
    numer = np.matmul(left_product_view, betas[:, :, np.newaxis])
    numer_view = numer.view().reshape(M, P)

    # Calculate result betas and standard errors
    new_betas = denom_recip_view * numer_view
    new_beta_ses = np.sqrt(denom_recip_view)

    return new_betas, new_beta_ses