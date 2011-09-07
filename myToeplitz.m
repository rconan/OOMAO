function t = toeplitz(c,r)
%TOEPLITZ Toeplitz matrix.
%   TOEPLITZ(C,R) is a non-symmetric Toeplitz matrix having C as its
%   first column and R as its first row.   
%
%   TOEPLITZ(R) is a symmetric Toeplitz matrix for real R.
%   For a complex vector R with a real first element, T = toeplitz(r) 
%   returns the Hermitian Toeplitz matrix formed from R. When the 
%   first element of R is not real, the resulting matrix is Hermitian 
%   off the main diagonal, i.e., T_{i,j} = conj(T_{j,i}) for i ~= j.
%
%   Class support for inputs C,R:
%      float: double, single
%
%   See also HANKEL.

%   Revised 10-8-92, LS - code from A.K. Booer.
%   Copyright 1984-2008 The MathWorks, Inc. 
%   $Revision: 5.11.4.3 $  $Date: 2008/06/20 08:00:11 $

% if nargin < 2,
%   c(1) = conj(c(1)); r = c; c = conj(c); % set up for Hermitian Toeplitz
% else
%   if r(1) ~= c(1)
%     warning('MATLAB:toeplitz:DiagonalConflict',['First element of ' ...
%            'input column does not match first element of input row. ' ...
%            '\n         Column wins diagonal conflict.'])
%   end
% end
%
r = r(:);                               % force column structure
p = length(r);
m = length(c);
x = [r(p:-1:2) ; c(:)];                 % build vector of user data
cidx = (0:m-1)';
ridx = p:-1:1;
t = cidx(:,ones(p,1)) + ridx(ones(m,1),:);  % Toeplitz subscripts
t = x(t);                                   % actual data

