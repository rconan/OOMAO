function H = bilinearSplineInterp(xo,yo,do,xi,yi)
%% BILINEARSPLINEINTERP Bilinear interpolation matrix
%
% H = bilinearSplineInterp(xo,yo,do,xi,yi,d); computes the sparse matrix H
% to perform the bilinear interpolation zi = H*zo, where zo and zi are
% defined on the meshes [xo;yo] and [xi;yi], respectively. do is the mesh
% step size of [xo;yo]

xo = xo(:)';
yo = yo(:)';
xi = xi(:);
yi = yi(:);

ni = length(xi);
if ni~=length(yi)
    error('xi and yi must have the same length.')
end
no = length(xo);
if no~=length(yo)
    error('xo and yo must have the same length.')
end

u = bsxfun(@minus,xi,xo)/do;
v = bsxfun(@minus,yi,yo)/do;

H = linearSpline(u).*linearSpline(v);

% [i,j,s] = find(H);
% [m,n]   = size(H);
% as      = abs(s);
% index   = as > eps(max(as)); % round-off error filter
% 
% H = sparse(i(index),j(index),s(index),m,n);

end

function y = linearSpline(x)
%% LINEARSPLINE Linear spline function
%
% y = linearSpline(x) computes the function y = 1 - |x| for |x|<1 and y = 0
% elsewhere

[m,n] = size(x);
x     = abs(x);
index = x < 1;
[i,j] = find(index);
s     = 1 - x(index);
y = sparse(i,j,s,m,n);

end
