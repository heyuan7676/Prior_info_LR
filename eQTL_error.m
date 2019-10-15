function [] = cross_validation()

% matrixeQTL error

for chr = 3:22
    [snps, genes, pairwise_features] = readin(chr);
    SNPs = snps;
    GENEs = genes;
    beta_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_beta.txt',chr);
    beta = importdata(beta_fn, '\t', 1);
    beta = beta.data;
    FDR_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_FDR.txt',chr);
    FDR = importdata(FDR_fn, '\t', 1);
    FDR = FDR.data;
    beta = beta .* (FDR<0.05);
    beta = zeros(size(beta));
    linear_error = norm(GENEs - beta * SNPs,'fro');
    fprintf('chr%d, error using betas =0 is %2.5f\n', chr, linear_error);
end

