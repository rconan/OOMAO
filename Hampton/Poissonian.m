function [dist, n] = Poissonian(lamda)
dist = ones(6*lamda,1);
dist(1) = exp(-lamda);
%for ki = 2:4*lamda
ki = 1;
k = 1;
while ki < lamda || dist(ki) > eps
    k = ki;
    ki = ki+1;
    dist(ki) = dist(k)*(lamda/k);
%         dist(ki) = (lamda^k)*(exp(-lamda))/factorial(k);
end
dist = dist(1:k);
n = 0:(size(dist,1)-1);
n = n(:);
    