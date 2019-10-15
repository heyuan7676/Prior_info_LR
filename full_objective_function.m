function obj = full_objective_function(W, Theta, S, pairwise_features, beta, C0, C1, lambda2, D, E, method)

    llk = logdet(Theta) - trace(S * Theta);
   
    [lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1); 
    if method == 3
        weight_penalty = sum(sum(lambda1 .* W.^2 ));  % l2 norm
    else
        weight_penalty = sum(sum(lambda1 .* abs(W)));  % l1 norm
    end

    % Theta = Theta;
    % Theta = ((eye(size(Theta)) - 1) .* -1) .* Theta;
    precision_penalty = sum(sum(lambda2 .* abs(Theta))); % l1 norm on off-diagonal elements

    beta_penalty = E * sum(beta.^2);
	
    obj = llk - precision_penalty - weight_penalty - beta_penalty;
