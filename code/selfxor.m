function [x,dbin] = selfxor(d,N,g)

if nargin < 2
    g = 1;
end

dbin = de2bi(d,N);
% dbin = dbin(end:-1:1); % change from little to big endian


% make sure can be divisible
if (mod(N,g))
    error('Cannot be grouped');
end

dbin = reshape(dbin,N/g,g);
% grps = N/g;
x = zeros(1,g);
for i=1:g
    if (mod(sum(dbin(:,i)),2))
        x(i) = 1;
    end
end