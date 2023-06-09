X = [1, 1; 2, 3; 4, 2; 5, 5];



n = size(X,1);
Xbar = X - ones(size(X)) * diag(mean(X));
Xcov = Xbar' * Xbar/n;
A = inv(Xcov);

V = X(2:n,:) - X(1:n-1,:);
t = diag(V * A * V') .^ (1/2);

V_2 = X(3:n,:) - X(1:n-2, :);
t_2 = diag(V_2 * A * V_2');

theta = zeros(n-1, 1);

for j = 2:n-1
    theta(j) = min(pi - acos( (t(j-1)-2 + t(j)-2 - ...
t_2(j-1))/(2*t(j)*t(j-1)) ) , pi/2); 
end

h = zeros(n - 1, 1);
disp(h);
h(1) = t(1)*( 1 + (1.5*theta(2)*t(2))/(t(1) + t(2)) );
disp(h(1));
for j = 2:n-2
h(j) = t(j) * ( 1 + (1.5*theta(j)*t(j-1))/(t(j-1)+t(j)) + ...
(1.5*theta(j+1)*t(j+1))/(t(j)+ t(j+1)) );
end 
disp(h(j));
h(n-1) = t(n-1) * ( 1 + (1.5*theta(n-1)*t(n-2))/(t(n-2)+t(n-1)) );

disp(h(n - 1));
disp(h);
h = [0;h];
disp(h);
h = cumsum(h);
disp(h);
h = h / h(n);
disp(h);
%end aff_angle
