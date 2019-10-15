function[W, beta, lambda1, Theta, all_objs] = callGLasso(chr, permute)

    chr = str2num(chr);
    permute = str2num(permute);

    % genes: genes x samples
    [snps, genes, pairwise_features] = readin(chr);
    [n_genes, n_samples]= size(genes);

    if permute
        %%% permute the samples
        fprintf('Permuted result\n');
        idx = randperm(n_samples);
        genes  = genes(:,idx);
    end

    y = genes';
    EmpricalS = cov(y);

    lambda2list = [0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5,1.0];
    for lambda2 = lambda2list
        [Theta, Sigma, opt, time, it, dGap] = QUIC('default', EmpricalS, lambda2, 1e-6, 1,  100);
        % calculate llk
        fX = log(det(Theta)) - trace(Theta * EmpricalS) - norm(lambda2.*Theta, 1);
        fprintf('chr%i, The llk is %2.5f for lambda2 = %2.5f, result in %2.5f non-zero off-diagonal entries in Theta\n', chr, fX, lambda2, (sum(sum(Theta ~=0)) - n_genes) /(n_genes * n_genes) );
    end
