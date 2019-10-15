function[W, beta, Theta, all_objs, testErr, Rsq2] = cross_sCGG(snps, genes, pairwise_features, C0, C1, E, lambda2, method, max_iter, cross_vali, idx, k, output_log_path, save_out, out_fn)
    % snps: SNPs x samples
    % genes: genes x samples
    % pairwise_features: cell array n_features x genes x SNPs
    diary(output_log_path);
    
    tol = 0.1;

    assert(size(snps, 2) == size(genes, 2));
    assert((sum(sum(isnan(snps))) == 0) & (sum(sum(isnan(genes))) == 0));
    n_pairwise_features = length(pairwise_features);
    for i = 1:n_pairwise_features
		assert(sum(sum(isnan(pairwise_features{i}))) == 0);
    end

    [n_snps, n_samples] = size(snps);
    n_genes = size(genes, 1);

    assert((size(pairwise_features{1}, 2) == n_snps) && (size(pairwise_features{1}, 1) == n_genes))
    not_converged = 1;


    % cross validation or not
    if cross_vali
        save_snps = snps;
        save_genes = genes;
        snps = snps(:,idx~=k);
        genes = genes(:,idx~=k);
    end


    % derive the training data

    x = snps';    % n x q
    y = genes';
    Cx = cov(x);   % q x q 
    Cy = cov(y);
    Cxy = bsxfun(@minus,x,mean(x))'*bsxfun(@minus,y,mean(y))/(size(x,1)-1);
    EmpricalS = Cy;

    % for GLasso, penalty only on the off-diagnol elements
    % lambda2=lambda2*(1-eye(n_genes));

    % initialize beta
    beta = ones(1, n_pairwise_features);
    beta = [-10,1];
    [lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1);
    
    % initialize W
    W = zeros(n_genes, n_snps);
    
    % initialize Theta
    SW = Cy- W * Cxy - Cxy' * W' + W * Cx * W';  % residual covariance 
    assert (sum(SW(:) ~= Cy(:)) == 0);
    Theta = diag(1./diag(Cy));
    Sigma = diag(diag(Cy));
    % [Theta, Sigma, opt, time, it, dGap] = QUIC('default', SW, lambda2, 1e-6, 0, 100);
    
 
    all_objs = [];
    iter = 0;
    D = 0;

    obj1 = full_objective_function(W, Theta, SW, pairwise_features, beta, C0, C1,lambda2, D, E , method);
    % fprintf('At iteration %i, objective is %2.5f\n', iter, obj1)

    while(not_converged)
	iter = iter + 1;

	% update Theta
        % fprintf('Run Graphical Lasso\n');
        [Theta, Sigma, opt, time, it, dGap] = QUIC('default', SW, lambda2, 1e-6, 0, 100, Theta, Sigma);
        Rsq2 = scale_free_toplogy(Theta);
        obj2 = full_objective_function(W, Theta, SW, pairwise_features, beta, C0, C1, lambda2, D, E, method);
        % fprintf('There are %d non-zero off-diagonal entries in Theta (%2.3f)\n', sum(sum(Theta~=0))-n_genes, (sum(sum(Theta~=0))-n_genes) / (n_genes * n_genes));
        % fprintf('Scale free toplogy R2 = %2.3f\n', Rsq2);

        % update W/Gamma and SW
        % fprintf('Update W\n');
        W = compute_snp_gene_weights(snps, genes, lambda1, 0, W, Cx, Cxy, Theta, method);
        SW = Cy- W * Cxy - Cxy' * W' + W * Cx * W';  % re-compute the residual covariance
        if cross_vali
            test_snps = save_snps(:,idx==k);
            test_genes = save_genes(:,idx==k);
        else
            test_snps = snps;
            test_genes = genes;
        end
        pre_y = W * test_snps;
        testErr = norm(pre_y - test_genes, 'fro');
        obj3 = full_objective_function(W, Theta, SW, pairwise_features, beta, C0, C1,lambda2, D, E, method);
        % fprintf('There are %d non-zero entries in W (%2.3f)\n', sum(sum(W~=0)), sum(sum(W~=0)) / (n_snps * n_genes));
        % fprintf('Frobenious norm of residuals: %2.3f\n', testErr);


        % update beta
        % fprintf('Update beta\n');
        beta = compute_feature_weights(pairwise_features, W, C0, C1, E, beta, method);
        [lambda1, p_regulator_causal] = compute_snp_priors(pairwise_features, beta, C0, C1);
        obj4 = full_objective_function(W, Theta, SW, pairwise_features, beta, C0, C1, lambda2, D, E, method);
        % fprintf('Mean value of lambda1: %2.3f\n', mean(lambda1(:)));

        
	% fprintf('Difference between obj2 and obj1, obj3 and obj2, obj4 and obj3: %2.3f, %2.3f, %2.3f\n', obj2-obj1, obj3-obj2, obj4 - obj3);

	all_objs(end + 1) = obj4;
	% fprintf('At iteration %i, objective is %2.3f\n', iter, obj4)


        % stop early for non-approatiate panelty setting
        if (sum(Theta(:)~=0) < 20 + n_genes)  | (sum(W(:)~=0) < 20)  | (abs(obj4) > 1e40)
            testErr = 1e40;
            break
        end

        % stop when stopping criteria met
	if (abs(obj1 - obj4) < tol)
		fprintf('Objective converged at iteration %i; quitting.\n', iter);
		break
	end
	if (iter > max_iter)
		fprintf('Maximum iterations (%i) reached\n', max_iter);
		break
	end
        obj1 = obj4;
    end


    if (testErr > 1e20)  
        if sum(Theta(:)~=0) < 20 + n_genes
            fprintf('Model not converged: Too few non-zero elements in Theta\n');
        elseif sum(W(:)~=0) < 20
            fprintf('Model not converged: Too few non-zero elements in Gamma\n');
        end
    elseif cross_vali == 0
        fprintf('There are %d non-zero entries (%2.5f), %d eGenes in W (%2.5f), error=%2.5f\n', sum(sum(W~=0)), sum(sum(W~=0)) / (n_snps * n_genes), sum(sum(W~=0,2) ~=0), sum(sum(W~=0,2) ~=0)/n_genes, testErr);
        fprintf('There are %d non-zero off-diagonal entries in Theta (%2.5f), Rsq2 = %2.5f\n', sum(sum(Theta~=0))-n_genes, (sum(sum(Theta~=0))-n_genes) / (n_genes * n_genes), Rsq2);
    end

    if save_out
        save(out_fn,'W','beta','lambda1','Theta','all_objs');
    end

end



