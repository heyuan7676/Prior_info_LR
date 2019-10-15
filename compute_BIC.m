function BIC = compute_BIC(W, Theta, S)
    llk = log(det(Theta)) - trace(S * Theta);
    n_samples = size(W,2);
    sn = sum(sum(Theta~=0));
    kn = sum(sum(W~=0));
    BIC = -2 * llk + (sn/2 + kn + 2) * log(n_samples);
