function [f, grad] = betaObjectiveFunction(beta, pairwise_features, W, C0, C1, E, method)
        [lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1);
        % log_term = sum(sum(-log(lambda1)));
        log_term = 0;
        if method == 3
            norm_term = sum(sum((W.^2) .* lambda1));
        else
            norm_term = sum(sum(abs(W) .* lambda1));
        end
        regularization_term = E * sum(beta.^2);
        f = log_term + norm_term + regularization_term;

        grad = zeros(size(beta));
        for k = 1:length(grad)
                dP_dbetak = p_regulator_causal .* (1 - p_regulator_causal) .* pairwise_features{k};
                dC_dbetak = C1 * dP_dbetak - C0 * dP_dbetak;
                % grad(k) = 2 * E * beta(k) - sum(sum(dC_dbetak ./ lambda1)) + sum(sum(abs(W) .* dC_dbetak));
                grad(k) = 2 * E * beta(k) + sum(sum(abs(W) .* dC_dbetak));
        end

end
