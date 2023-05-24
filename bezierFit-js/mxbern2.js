/* 
    mxbern2.js

    function [B] = mxbern(2, d)

    // This function creates a bernstein matrix of degreee d using 
    // the values in the column vector t. 

    // berstein matrix is a generalized Vandermonde matrix whose (i,j) entry is 
    // B_{j-1}^d(t_j). The vector t must be a column vector with values between
    // 0 and 1

    [n m] = size(t);
    if (m ~= 1)
        error('t must be a column vector.');
    end

    // check if nodes are withing [0, 1]
    
    if min(t) < 0 | max(t) > 1
        error('nodes are not within [0, 1]')
    end

    // build bernstein matrix
    ct = 1 - t;
    B = zeros(n, d+1);

    for i = 0 : d
        B (:, i+1) = (t.^i).*(ct.^(d-1));
    end

    // if d > 22 we form matrix differently to avoid roundoff. 

    if d < 23
        B = B*diag( [1 cumprod(d: -1:1)./cumprod(1:d)]);
    else
        B = B*diag(diag(fliplr(pascal(d+1))));
    end
*/