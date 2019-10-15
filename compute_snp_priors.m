function [lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1)
        %This computes the C_ms (n_snps x n_modules). Larger C -> smaller SNP weights. Checked. 
        n_pairwise_features = length(beta);
        summed_weighted_features = zeros(size(pairwise_features{1}));
        for i = 1:n_pairwise_features
                summed_weighted_features = summed_weighted_features + beta(i) * pairwise_features{i};
        end
        p_regulator_causal = sigmoid(summed_weighted_features);
        lambda1 = C1 * p_regulator_causal + C0 * (1 - p_regulator_causal);
end

function x = sigmoid(z)
        x = (1 + exp(-z)) .^ (-1);
end
