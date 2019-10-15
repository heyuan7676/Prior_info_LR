function [] = callsCGG(chr)

chr = str2num(chr);
snp_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_genotypes.txt', chr);
snps = importdata(snp_fn, '\t', 1);
samples = snps.textdata(1,:);


distance_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_DISTANCE.txt',chr);
distance = importdata(distance_fn, '\t', 1);

contacting_fn = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_Feature_contacting_value.txt',chr);
contacting = importdata(contacting_fn, '\t', 1);


[n_snps, n_samples] = size(snps.data);
[n_genes, n_snps] = size(contacting.data);

% make sure that the samples, snps and genes align

assert(n_snps - sum(strcmp(snps.textdata(:,1)',distance.textdata(1,:))) <=1 );

% read in 

snps = snps.data;
n_pairwise_features = 2;

pairwise_features = cell(n_pairwise_features);
pairwise_features{1} = distance.data; % distance
pairwise_features{2} = contacting.data; % contacting values


% normalize

for k = 1:n_pairwise_features
    pairwise_features{k} = (pairwise_features{k} - min(min(pairwise_features{k})) ) / (max(max(pairwise_features{k})) - min(min(pairwise_features{k})));
end

snps = snps - repmat(mean(snps,2), 1, size(snps,2));


% initialize the parameters

E = 0.001;
C1 = 0.5;
C0 = 5;
beta = [-0.2, 0.25];

% generate the penalty weights
[lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1);

% generate weights matrix
Gamma = normrnd(0, lambda1);
Gamma(Gamma < quantile(Gamma(:), 0.9) ) = 0;
mu = Gamma * snps;

% generate sparse concentration matrix
sigma = sprandsym(n_genes, 0.3, rand(1,n_genes));
sigma(logical(eye(n_genes))) = 1;

% generate Y matrix
Y = mvnrnd(mu',inv(sigma));
Y = Y'; 



% save

outfile = sprintf('/scratch1/battle-fs1/heyuan/Predict_target_gene/GTEx/for_Glasso/GTEx_chr%d_pairs_activeSNP_genes_corrected.txt', chr); 
fid  = fopen(outfile, 'w');
for ii = 1:size(Y,1)
    fprintf(fid,'%g\t',Y(ii,:));
    fprintf(fid,'\n');
end
fclose(fid);


