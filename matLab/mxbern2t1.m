t = [
    0, 0.00781250000000000, 0.0156250000000000, 0.0234375000000000, 0.0312500000000000, ...
    0.0390625000000000, 0.0468750000000000, 0.0546875000000000, 0.0625000000000000, ...
    0.0703125000000000, 0.0781250000000000, 0.0859375000000000, 0.0937500000000000, ...
    0.101562500000000, 0.109375000000000, 0.117187500000000, 0.125000000000000, ...
    0.132812500000000, 0.140625000000000, 0.148437500000000, 0.156250000000000, ...
    0.164062500000000, 0.171875000000000, 0.179687500000000, 0.187500000000000, ...
    0.195312500000000, 0.203125000000000, 0.210937500000000, 0.218750000000000, ...
    0.226562500000000, 0.234375000000000, 0.242187500000000, 0.250000000000000, ...
    0.257812500000000, 0.265625000000000, 0.273437500000000, 0.281250000000000, ...
    0.289062500000000, 0.296875000000000, 0.304687500000000, 0.312500000000000, ...
    0.320312500000000, 0.328125000000000, 0.335937500000000, 0.343750000000000, ...
    0.351562500000000, 0.359375000000000, 0.367187500000000, 0.375000000000000, ...
    0.382812500000000, 0.390625000000000, 0.398437500000000, 0.406250000000000, ...
    0.414062500000000, 0.421875000000000, 0.429687500000000, 0.437500000000000, ...
    0.445312500000000, 0.453125000000000, 0.460937500000000, 0.468750000000000, ...
    0.476562500000000, 0.484375000000000, 0.492187500000000, 0.500000000000000, ...
    0.507812500000000, 0.515625000000000, 0.523437500000000, 0.531250000000000, ...
    0.539062500000000, 0.546875000000000, 0.554687500000000, 0.562500000000000, ...
    0.570312500000000, 0.578125000000000, 0.585937500000000, 0.593750000000000, ...
    0.601562500000000, 0.609375000000000, 0.617187500000000, 0.625000000000000, ...
    0.632812500000000, 0.640625000000000, 0.648437500000000, 0.656250000000000, ...
    0.664062500000000, 0.671875000000000, 0.679687500000000, 0.687500000000000, ...
    0.695312500000000, 0.703125000000000, 0.710937500000000, 0.718750000000000, ...
    0.726562500000000, 0.734375000000000, 0.742187500000000, 0.750000000000000, ...
    0.757812500000000, 0.765625000000000, 0.773437500000000, 0.781250000000000, ...
    0.789062500000000, 0.796875000000000, 0.804687500000000, 0.812500000000000, ...
    0.820312500000000, 0.828125000000000, 0.835937500000000, 0.843750000000000, ...
    0.851562500000000, 0.859375000000000, 0.867187500000000, 0.875000000000000, ...
    0.882812500000000, 0.890625000000000, 0.898437500000000, 0.906250000000000, ...
    0.914062500000000, 0.921875000000000, 0.929687500000000, 0.937500000000000, ...
    0.945312500000000, 0.953125000000000, 0.960937500000000, 0.968750000000000, ...
    0.976562500000000, 0.984375000000000, 0.992187500000000, 1
];
t = t(:);
d = 3;
[n, m] = size(t);
ct = 1 - t;
B = zeros(n,d+1);

for i = 0:d
B(:,i+1)= (t.^i).*(ct.^(d-i));
end
disp(d);
j = [1 cumprod(d:-1:1)./cumprod(1:d)];
if d < 23
B = B*diag( [1 cumprod(d:-1:1)./cumprod(1:d)] );
else
B = B*diag(diag(fliplr(pascal(d+l))));
end
disp(B);