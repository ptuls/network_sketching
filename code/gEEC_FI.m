function [FI,CRB] = gEEC_FI(theta,b,li,K)
%% Comparing Fisher information for gEEC scheme for sampling without
%  replacement versus sampling with replacement approximation. Here,
%  immunity was assumed. The Fisher information is computed for any single 
%  level with immunity. 
%
%  Input
%  theta: channel bit flip probability
%  b    : number of bits of subsketch
%  li   : number of levels of subsketches
%  K    : alphabet size
%
%  Output
%  FI   : total Fisher information of all subsketches
%  CRB  : Cramer-Rao lower bound of all subsketches combined
%
%  Paul Tune Mon 04 Aug 2014
%  paul.tune@adelaide.edu.au
%
%  Last updated: Mon 04 Aug 2014
%  

if nargin < 4
    K = 2;
end

if theta <= 0 || theta >= 1
    error('theta must be between 0 and 1');
end

if b < 0
    error('Number of bits of subsketch must be positive');
end

if li < 1
    error('Number of levels of subsketch must at least 1');
end

if K < 2
    error('Alphabet size must be at least 2');
end

% Modified version: compute via circular convolution
% much faster
m = 1-theta;
md = -1;
if (K > 1)
    m = [m theta/2];
    md = [md 0.5];
    if (K > 2)
        m = [m zeros(1,K-3) theta/2];
        md = [md zeros(1,K-3) 0.5];
    end
end

% Construct probabilities and differential
f = zeros(1,K);
fd = f;
f(1) = 1;

FI = 0;
for j=1:li
    % watch out for order here
    fd = cconv(fd,m,K) + cconv(f,md,K);
    % each time in the recursion, am computing level j
    f = cconv(f,m,K);
    % clean up noise from convolution
    f(f < 0) = 0;

    % prevent division by 0
    % sum all contributions from all levels
    FI = FI+sum(fd(f > 0).^2./f(f > 0));
end

% total based on number of bits (assuming all levels have same size)
FI = b*FI;
CRB = 1/FI;