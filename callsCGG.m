function [] = callsCGG(chrlist, permute, seq)

chrlist = str2num(chrlist);
permute = str2num(permute);

lambda2 = 0.0550;
C1 = 0.1100;
C0 = C1 * 5;


fprintf('Perform sCGG for chromosome: %g\n', chrlist);

chr = chrlist(1);
[snps, genes, pairwise_features] = readin(chr);
SNPs = snps;
GENEs = genes;
PF = pairwise_features;

if length(chrlist) > 1
    for chr = chrlist(2:end)
        [snps, genes, pairwise_features] = readin(chr);
        SNPs = [SNPs; snps];
        GENEs = [GENEs; genes];
        for k = 1: length(PF)
            PF{k} = [PF{k} zeros(size(PF{k},1), size(snps,1)); zeros(size(genes,1), size(PF{k},2)) pairwise_features{1}];
        end
    end
end

[n_snps, n_samples] = size(SNPs);
n_genes = size(GENEs, 1);
fprintf('Running sCGG with %i snps, %i genes, %i samples\n', n_snps, n_genes, n_samples)


if permute
    %%% permute the samples
    fprintf('Permuted result for chr%i\n', chr);
    permuted_sample_idx = randperm(n_samples);
    SNPs = SNPs(:,permuted_sample_idx);
    output_log_path=sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/run_callsCGG_output_chr%d_permuted.txt',chr);
    out_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_results_permuted_%d.mat', chr, seq);
else
    SNPs = SNPs;
    output_log_path=sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/run_callsCGG_output_chr%d.txt',chr);
    out_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_results.mat', chr);
end


savefn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/train_test_idx/permute_idx.mat');
load(savefn)


E = 0.001;
max_iter=100;

tic
fprintf('Update W using Rothmans method\n');
[W, beta, Theta, all_objs, testErr, Rsquare] = sCGG(SNPs, GENEs, PF, C0, C1, E, lambda2, 2, max_iter, 0, idx, 0, output_log_path, 1, out_fn);
fprintf('chr%d (%d snps, %d genes), C1=%2.4f, C0=%2.2f, lambda2=%2.4f:  prediction error =%2.5f, Rsquare=%2.5f\n', chr, n_snps, n_genes, C1, C0, lambda2, testErr, Rsquare);
fprintf('\n');
toc

beta


