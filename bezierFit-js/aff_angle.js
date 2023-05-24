/* 
    aff_angle.js

    function [h] = aff_angle(X)

    this function retruns affine invariant vector of nodes for a given set of ordered data X.
    it is assumed that the user sends this function the data arranged in a 
    (n x 2) matrix ordered from row one to row n

    obtaining the covariance matrix A is from affknt.m 
    copyrighted 1994 Carlos F. Borges
    the affine invariant angle method of obtaining nodes -
    Thomas A. Foley and Gregory M. Nielson "Knot Selection for Parametric Spline Interpolation"
    in "Mathematical methods in computer aided geometric design"

    n = size(X, 1);
    Xbar = X - ones(size(X)) * diag(mean(X));
    Xcov = Xbar' * Xbar/n;
    A = inv(Xcov);

    obtain node spacing values using metric
        t_i = M[X] (X_i, X_(i+1))

    V = X(2:n, :) - X(1:n-1, :);

    t = diag(V * A * V') .^(1/2);

    obtain values for
        M^2[X](X_(i-1), X(i+1))
    
    V_2 = X(3:n, :) - X(1:n-2, :);
    t_2 = diag(V_2 * A * V_2');

    get theta_i values

    theta = zeros(n-1, 1);

    for j = 2: n-1
        theta(j) = min (pi - acos((t(j-1) ^2 + t(j)^2 ~ ...
                    t_2(j-1))/(2*t(j)*t(j - 1)) ) , pi/2); 
    end

    // Obtain the affine invariant angle node spacing values h_i

    h = zeros(n - 1, 1);

    h(1) = t(1) * (1 + (1.5*theta(2) * t(2)) /(t(1) + t(2)) );

    end

    h(n -1) = t(n-1) * (1 + (1.5 * theta(n - 1) * t (n -2)) / (t (n -2) + t(n -1)));

    // the nodes are now spacing values, so normalize them so they are within [0, 1],
    // with first data point being associated with value zero and teh last data point with
    // the value one

    h = [0 ; h];

    h = cumsum(h);

    h = h / h(n);


    


*/