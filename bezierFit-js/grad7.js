// function [p, t, info] = grad7(d, deg, stop)

/* 
    Input Args:
        d: matrix of data points
        deg : degree of curve the user wants to fit
        stop: stopping criteria when 10^(stop)
    
    Output Args:
        p: control point for best fit Bezier
        t: nodes for best fit Bezier
        info: 2x1 vector which has final norm of the residual and the number of iterations to convergence

        
    step 1. 
    if (ishold)
        hold_was_off = 0;
    else
        hold_was_off = 1;
    end;

    // i: number of data points j: number of control points
    // aff_angle(d) returns initial set of nodes (a vector of i elements)
    i = size(d, 1);
    j = deg + 1;
    t = aff_angle(d);

    // bez_mat: (i x j) matrix, determined by nodes and degree of curve
    // j x 2 matrix of control points 'p'. aka. linear least squares problem
    // 'p' determined by p = pinv(bez_mat) * d

    bez_mat = mxbern2(t, deg);
    p = bez_mat \ d; // == p = pinv(bez_mat) * d

    // with control points p, we can find curve and points on the curve associated with initial vector t.
    // y is the i x 2 matrix of points on Bezier curve associated with inital vector t
    // t1 is a closely spaced set of parameter values which produce bezier matrix
    // bez_mat_1 will give enough points on fitted bezier curve for matlab to plot a smooth looking curve
    // y1 is matrix of these poitns

    y = bez_mat * p;
    t1 = [0:1 / 128:1]';
    bez_mat_1 = mxbern2(t1, deg); 
    y1 = bez_mat_1 * p;

    figure
    plot(y1 (: , 1), y1(:, 2))
    hold on
    plot(p (: , 1), p(: , 2), '*')
    plot(y (: , 1), y(: , 2), 'o')
    plot(d (: , 1), d(: , 2), '+')
    axis('equal')
    legend('-', 'Fitted Curve', '*', "Control Points',...
            'o', 'Initial Times', '+', 'Data Points')
    pause

    step 2.

    iter = 0;
    resid_old = 0;
    resid_new = bez_mat * p - d;
    tic

    // B(t) * p = d
    //each iterationg of the loop produces a new vector t, by solving nonlinear least squares problem
    
    while norm(resid_new - resid_old)/ max(1, norm(resid_new)) > 10 ^ (stop)
    
    // where 'p' is the vector of control points
    // 'd' the matrix of data points
    // B(t) is bernstein matrix with unkown parameter values
    // the method used in this function is the Gauss-Newton method
    // the residual is R(t) = B(t) * p - d and the Jacobian matrix 'J'
    // be such that J-i,j = dB(t_i).dt_j
    // AKA  (i, j)^th element of J is the slope along the Bezier curve at the i^th parmeter value
    // with respect to the j^th parameter value
    // 
    // the Gauss-Newton method says taht the change in paramter values which will minimize the residual given
    // by delta_t = -inv(J' * J) * J; * R
    // where 'J' and 'R' are evaluated at the current parameter values

    
    // computes the gradient at a point on the bezier curve
    // this slope (deriv) is found by multiplying Bernstein matrix for a degree - 1 curve
    // by the forward difference matrix of control points, then multiplying by the degree of the original bezier
    deriv = deg * mxbern2(t, deg-1) * (p(2:j, ) - p(1:j-1,:));

    // y-values of the residual and append them to bottom of x-value (resid(:))
    // matrix 'J' is x or y value of slope at each parameter's point
    // most entries for 'J' are zero

    // 't' are the new nodes given by 't = t - delta_t'
    // we can form 'delta_t' using less computer time and flops
    t = t - (deriv(:, 1). *resid_new(:, 1) + deriv(:, 2). *resid_new(:, 2))...
            ./ [deriv(:,1).^2 + deriv(:,2).^2];
    

    // with new vector t, check to make sure values are between 0 and 1
    // following rescaling of the nodes also results in the endpoints being associated with
    // nodes t_1=0 and t_m=1

    t = -min(t)*ones(i, 1) + t;
    t = t/max(t);

    // With ordered nodes we now get new control points 
    // to reproduce the points 'tau' on the curve for the next iteration of the while loop
    // if the condition is met we can use bez_mat matrix for the final plot
    // a new value for resid_new is computed and the iteration counter is updated

    bez_mat = mxbern2(t, deg);
    P = bez_mat \ d; 

    resid_old = resid_new;
    resid_new = bez_mat * p - d;

    iter = iter + 1;

    end
    toc

    // we now have the best fit vector t and associated matrix p of control points
    // y are the points on bezier curve associated with vector t
    // y1 are the closely space points on the bezier curve

    // y = bez_mat * p
    // y1 = bez_mat_1 * p;

    figure
    plot(y1(:,1),y1(:,2))
    hold on
    plot(p(:,1),p(:,2),'*')
    plot(y(:,1),y(:,1),'o')
    plot(d(:, 1),d(:,2), '+')
    axis('equal')

    legend('-', 'Fitted curve', '*', 'Control Points', ...
            'o', 'Nodes' , '+', 'Data Points')

    if (hold_was_off)
        hold off;
    end

    info = [norm(bez_mat*p-d), iter];
*/


