function [] = BIC_validation(chr, method, permute)

chr = str2num(chr);
method = str2num(method);
permute = str2num(permute);

[snps, genes, pairwise_features] = readin(chr);
SNPs = snps;
GENEs = genes;
PF = pairwise_features;
[n_snps, n_samples] = size(snps);
n_genes = size(genes, 1);
fprintf('Compute BIC using sCGG for chr%i, with %i snps, %i genes, %i samples\n', chr, n_snps, n_genes, n_samples);

if method == 3
    fprintf('Update W using l2 norm\n');
elseif method == 2
    fprintf('Update W using Rothmans method\n');
end


if permute
    %%% permute the samples
    fprintf('Permuted result for chr%i\n', chr);
    idx = randperm(n_samples);
    SNPs = SNPs(:,idx);
else
    SNPs = SNPs;
end



E = 0.001;
max_iter=100;

output_log_path=sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/run_callsCGG_output_chr%d.txt',chr);
out_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_results.mat', chr);



for lambda2 = [0.2, 0.5, 0.75]
    for C1 = [0.001, 0.005, 0.01]
        for C0 = [C1.*5, C1.*10]
            tic
            [W, beta, lambda1, Theta, all_objs, BIC] = BIC_sCGG(SNPs, GENEs, PF, C0, C1, E, lambda2, method, max_iter, output_log_path, out_fn);
            fprintf('chr%d (%d snps, %d genes), C1=%2.2f, C0=%2.2f, lambda2=%2.2f:  test error =%2.5f\n', chr, n_snps, n_genes, C1, C0, lambda2, errorSum);
            toc
            fprintf('\n');
        end
    end
end

