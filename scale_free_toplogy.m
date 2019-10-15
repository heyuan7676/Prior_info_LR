function Rsq2 = scale_free_toplogy(Theta)
        n_genes = size(Theta,2);
        degree = sum(Theta~=0);
        d = unique(degree);

        if sum(Theta(:)~=0) < 20 + n_genes | length(d) == 1
                Rsq2 = 0;
        else
                pd = zeros(1,length(d));
                for kd = 1:length(d)
                    pd(kd) = sum(degree == d(kd))/n_genes;
                end
                y = log(pd)';
                x = [ones(length(d),1) log(d)'];
                b = x \ y;
                ypred = x * b;
                Rsq2 = 1 - sum((y - ypred).^2)/sum((y - mean(y)).^2);
        end

end
