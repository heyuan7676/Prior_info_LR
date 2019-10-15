import pandas as pd
import numpy as np
import os

chromosome = int(os.environ['chr'])
beta = [-0.2,0.25]
c1 = 0.05
c2 = 5

PAIR_DIR = '/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/'

genes_chr = pd.read_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_genes_corrected.txt' % chromosome), sep='\t', index_col=0)
snps_chr = pd.read_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_genotypes.txt' % chromosome), sep='\t', index_col=0)

f = 'DISTANCE'
DIS = pd.read_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_Feature_%s.txt' % (chromosome, f)), sep='\t', index_col=0)

f = 'contacting_value'
contacting = pd.read_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_Feature_%s.txt' % (chromosome, f)), sep='\t', index_col=0)



### pick a random set of SNPs and Genes

SNPs = np.random.choice(list(snps_chr.index), size = len(snps_chr), replace=False)
Genes = np.random.choice(list(genes_chr.index), size = len(genes_chr), replace=False)


### get the X and F matrices

chromosome = int(os.environ['chr_out'])
X = snps_chr.loc[SNPs]
X.to_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_genotypes.txt' % chromosome), sep='\t')

F1 = DIS.loc[Genes][SNPs]
f = 'DISTANCE'
F1.to_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_Feature_%s.txt' % (chromosome, f)), sep='\t')

F2 = contacting.loc[Genes][SNPs]
f = 'contacting_value'
F2.to_csv(os.path.join(PAIR_DIR, 'GTEx_chr%i_pairs_activeSNP_Feature_%s.txt' % (chromosome, f)), sep='\t')

