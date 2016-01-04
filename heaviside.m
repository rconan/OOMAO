function out = heaviside( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

out = zeros(size(x));
out(x==0) = 0.5;
out(x>0) = 1;

end

