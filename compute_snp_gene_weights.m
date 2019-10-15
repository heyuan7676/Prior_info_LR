function W = compute_snp_gene_weights(snps, genes, lambda1, initialize, W, Cx,Cxy, Theta, method)
    % Optimize W
    [n_snps, n] = size(snps);
    n_genes = size(genes,1);
    
    % initialize W using separate lasso regressions
    if initialize
        W = zeros(n_genes, n_snps);
        for k = 1:n_genes
            A = snps';
            b = genes(k,:)'; 
            lambda = lambda1(1,1);
            rho = 0.0001;
            alpha = 1.0;
            [z, history] = lasso(A, b, lambda, rho, alpha);
            W(k,:) = z;
        end
    % update W
    else
        if method == 1
            ej = ones(n_snps,1);
            ei = ones(n_genes,1);
            gVector = 2 * (ej' * (Cxy * Theta) * ei + (ej'*Cx*ej) *(ei'*Theta*ei)*W - ej'*(Cx*W'*Theta)*ei);
            W = (abs(gVector - lambda1) > 0 ).*abs(gVector - lambda1).* sign(gVector) ./ (2*(ej'*Cx*ej) * (ei'*Theta*ei)); 
        elseif method == 2
            X = snps';
            Y = genes';
            lambda1 = lambda1';
            S = X' * X;
            H = X' * Y * Theta;
            W = W';
            %for r = 1:n_snps
            %    for c = 1:n_genes
            %        U(r,c) = S(r,:) * (W * Theta(:,c));
            %        A(r,c) = W(r,c) + (H(r,c) - U(r,c)) / (S(r,r) * Theta(c,c));
            %        B(r,c) = abs(A(r,c)) - n*lambda1(r,c) / (S(r,r) * Theta(c,c)) ;
            %        w_rc = sign(A(r,c)) * (B(r,c)>0) * B(r,c);
            %        W(r,c) = w_rc;
            %    end
            %end
            %W = W';
            
            % update W row by row (per snp), 45 times faster 
            % fails to do it column by column
            % probably because the effect on all the genes for one SNP is more similar than all SNPs' effect on one gene

            Denominator = diag(S) * diag(Theta)';
            for r = 1:n_snps
                Ur = S(r,:) * W * Theta;
                Ar = W(r,:) + (H(r,:)-Ur) ./ Denominator(r,:);
                Br = abs(Ar) - n.*lambda1(r,:) ./ Denominator(r,:);
                W(r,:) = sign(Ar) .* (Br>0) .* Br;
            end
            W = W';

        elseif method == 3
            % penalize l2 norm
            A = inv(Theta);
            B = snps * snps' ./ n;
            C = genes * snps' ./ n;
            result = pinv([A eye(n_genes)]) * C * pinv([eye(n_snps);B]);
            W = result(1:n_genes, 1:n_snps);
            W = W./lambda1;
       end
    end

