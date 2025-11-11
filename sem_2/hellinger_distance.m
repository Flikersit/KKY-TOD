%% === Hellinger distance function ===
function H2 = hellinger_distance(mu1, S1, mu2, S2)
    % Computes the squared Hellinger distance between two 2D Gaussian distributions
    % If the means are identical, the metric depends only on the covariance matrices
    
    logdet = @(A) 2 * sum(log(diag(chol(A, 'lower')))); % stable log(det)
    S = 0.5 * (S1 + S2);
    dmu = mu1 - mu2;
    
    % Hellinger^2 formula
    H2 = 1 - exp( 0.25 * (logdet(S1) + logdet(S2)) ...
                  - 0.5 * logdet(S) ...
                  - 0.125 * (dmu' * (S \ dmu)) );
end

