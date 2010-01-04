% Zernike spatial covariance

D = 10;
L0 = 30;
x0 = pi*D/L0;

fun = @(x,n1,n2,m) (x.^2+x0.^2).^(-11/6).*besselj(n1+1,x).*besselj(n2+1,x).*besselj(m,4.*x)./x;

tic
quadgk(@(x) fun(x,1,1,0), 0, Inf)
toc