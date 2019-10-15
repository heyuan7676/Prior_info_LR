function [snps, genes, pairwise_features] = readin(chr)

% chr = str2num(chr);
snp_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_genotypes.txt', chr);
snps = importdata(snp_fn, '\t', 1);

gene_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_genes_corrected.txt', chr);
genes = importdata(gene_fn, '\t',1);

distance_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_DISTANCE.txt',chr);
distance = importdata(distance_fn, '\t', 1);

contacting_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_contacting_value.txt',chr);
contacting = importdata(contacting_fn, '\t', 1);

% beta_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_beta.txt',chr);
% beta = importdata(beta_fn, '\t', 1);

[n_snps, n_samples] = size(snps.data);
n_genes = size(genes.data,1);



% make sure that the samples, snps and genes align

assert(n_samples - sum(strcmp(genes.textdata(1,:), snps.textdata(1,:))) <= 1);
assert(n_snps - sum(strcmp(snps.textdata(:,1)',distance.textdata(1,:))) <=1 );
% assert(n_genes - sum(strcmp(genes.textdata(:,1), distance.textdata(:,1))) <=1 );

% read in 

snps = snps.data;
genes = genes.data;

pairwise_features = cell(2);
pairwise_features{1} = distance.data; % distance
pairwise_features{2} = contacting.data; % contacting values

% beta = beta.data;

% normalize 
normalizeFeatures = 1;
if normalizeFeatures
    for k = 1:2
        pairwise_features{k} = (pairwise_features{k} - min(min(pairwise_features{k})) ) / (max(max(pairwise_features{k})) - min(min(pairwise_features{k})));
    end
end

snps = snps - repmat(mean(snps,2), 1, size(snps,2));
genes = genes - repmat(mean(genes,2), 1, size(genes,2));
