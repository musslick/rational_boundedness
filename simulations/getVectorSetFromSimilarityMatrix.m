function X = getVectorSetFromSimilarityMatrix(Sigma, nDim)

X = randn(nDim, size(Sigma,1));
X = bsxfun(@minus, X, mean(X));
X = X * inv(chol(cov(X)));
X = X * chol(Sigma);

end