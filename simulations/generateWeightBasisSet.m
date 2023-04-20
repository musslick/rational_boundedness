function weights = generateWeightBasisSet(nHiddenUnits, nTaskUnits, similarity)

n = nHiddenUnits;
k = nTaskUnits;
d = sqrt(nTaskUnits);
Sigma = eye(nTaskUnits, nTaskUnits);

% set up Sigma
for row =1:k
    for col = 1:k
        row_block = ceil(row/d);
        col_block = ceil(col/d);
        if(row_block == col_block)
            if(row ~= col)
                Sigma(row, col) = similarity;
            end
        else
            Sigma(row, col) = 0;
        end
    end
end

X = getVectorSetFromSimilarityMatrix(Sigma, n);

% norm each column vector
for col = 1:size(X,2)
    n = norm(X(:, col));
    X(:, col) = X(:, col) / n;
end

weights = X;

end
