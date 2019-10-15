function [] = cross_validation(chr, method)

chr = str2num(chr);
method = str2num(method);

[snps, genes, pairwise_features] = readin(chr);
SNPs = snps;
GENEs = genes;
PF = pairwise_features;
[n_snps, n_samples] = size(snps);
n_genes = size(genes, 1);
% fprintf('Running 5-fold cross validtion of sCGG for chr%i, with %i snps, %i genes, %i samples\n', chr, n_snps, n_genes, n_samples);

if method == 3
    fprintf('Update W using l2 norm\n');
elseif method == 2
    fprintf('Update W using Rothmans method\n');
end


idx = crossvalind('Kfold', 1:n_samples);



E = 0.001;
max_iter=100;

output_log_path=sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/run_callsCGG_output_chr%d.txt',chr);
out_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_results.mat', chr);


% matrixeQTL error
% beta_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_beta.txt',chr);
% beta = importdata(beta_fn, '\t', 1);
% beta = beta.data;
% linear_error = norm(GENEs - beta * SNPs,2);
% fprintf('chr%d, error from linear regression = %2.5f\n', chr, linear_error);


for lambda2 = [0.2, 0.5, 0.75]
    for C1 = [0.001, 0.005, 0.01]
        for C0 = [C1.*5, C1.*10]
            tic
            errorSum = 0;
            for k = 1:5
                [W, beta, lambda1, Theta, all_objs, testErr] = cross_sCGG(SNPs, GENEs, PF, idx, k, C0, C1, E, lambda2, method, max_iter, output_log_path, out_fn);
                errorSum = errorSum + testErr;
            end
            fprintf('chr%d (%d snps, %d genes), C1=%2.2f, C0=%2.2f, lambda2=%2.2f:  test error =%2.5f\n', chr, n_snps, n_genes, C1, C0, lambda2, errorSum);
            toc
            fprintf('\n');
        end
    end
end







