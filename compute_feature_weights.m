function beta = compute_feature_weights(pairwise_features, W, C0, C1, E, old_beta, method)
        options = optimoptions(@fminunc, 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient', true, 'Display', 'none');
        beta = fminunc(@(beta) betaObjectiveFunction(beta, pairwise_features, W, C0, C1, E, method), old_beta, options);
        %fprintf('Objective value of final beta')
        %[f, grad] = betaObjectiveFunction(beta, pairwise_features, W, C0, C1, E)
        %fprintf('Objective value of true beta')
        %[f, grad] = betaObjectiveFunction([10, 0, 0], pairwise_features, W, C0, C1, E)
end

