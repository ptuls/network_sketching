%% Comparing Fisher information content without replacement versus sampling 
%  with replacement approximation. Here, immunity was assumed. The Fisher 
%  information is computed for any single level with immunity. 
%
%  Paul Tune Mon 28 Jul 2014
%  paul.tune@adelaide.edu.au
%
%  Last updated: Mon 28 Jul 2014
%  
clear;

% bit flip probability i.e. theta. Range this from 'lo' to 'hi'
lo = 0.001;
hi = 0.6;
lim = 100;
theta = linspace(lo,hi,lim);

% NB: this is a two level model. To analyse sampling with replacement 
% vs sampling without replacement, we need to specify the number of
% non-zero bits in the packet. There are two parameters: n for the number
% of non-zero bits, N for the total number of bits in the packet we sample
% from.
%
% Model: packet contents -> sampled bits -> channel -> observations
%
% Compute Fisher information of received bits.

% Samples
b = 1;
% Number of non-zero bits
n = 300;
% Packet size in bits (the population in sampling terminology)
N = 1000;

FIwr_total = zeros(1,lim);
FIwor_total = zeros(1,lim);
comb = zeros(1,b+1);

%% Compute bit sampling probabilities
%  sampling with replacement 
p = n/N;
q = 1-p;
L = 2^b;

% sampling with replacement probabilities (for comparison)
p1 = [q p];

% sampling without replacement (need to find better way of calculating)
for j=0:b
    comb(j+1) = nchoosek(b,j);
end
h = hygepdf(0:b,N,n,b)./comb;

% case input is (0,0,...,0),(0,0,...,1),...,(1,1,...,1) in that order
p2 = h(sum(de2bi(0:(L-1)),2)+1);  % probability of all input sequences

tic
Td1 = [-1 1; 1 -1];
for i=1:lim
    t = [1-theta(i) theta(i)];

    %% Fisher information: sampling without replacement
    %  probability of output of the channel for every input sequence
    %  essentially a Bernoulli process
    T1 = [t; t(end:-1:1)];
    
    T2 = T1;
    for k=1:(b-1)
        T2 = kron(T2,T1);
    end   
    
    % calculating differential w.r.t. theta
    % derivative for every input sequence into the channel
    Td2 = zeros(L,L);
    for j=1:b
        if (j==1)
            P = Td1;
        else
            P = T1;
        end
        for k=2:b
            if (j==k)
                P = kron(P,Td1);
            else
                P = kron(P,T1);
            end
        end
        Td2 = Td2 + P;
    end
    
    f2 = diag(p2)*T2;
    fd2 = diag(p2)*Td2;
    
    id2 = find(f2 > 0);
    % compute Fisher information
    FIwor_total(i) = sum(fd2(id2).^2./f2(id2));
    
    %% Fisher information: sampling with replacement
    %  since i.i.d., so multiply with subsketch size 
    
    f1 = diag(p1)*T1;
    fd1 = diag(p1)*Td1;
    
    id1 = find(f1 > 0);
    % compute Fisher information
    FIwr_total(i) = b*(q-p)^2./((q*(1-theta(i)) + p*theta(i))*(q*theta(i) ...
    + p*(1-theta(i))));
end
toc

print_figures = 0;
save_data = 0;
fontsize = 18;
dgrey = [0.65 0.65 0.65];
linewidth = 2;
position = [0 0 14 12];

set(0, 'defaultTextInterpreter', 'latex'); 
figure(1)

h1 = semilogy(theta,theta.^2.*FIwor_total, 'color', dgrey, 'linewidth', linewidth);
hold on
h2 = semilogy(theta,theta.^2.*FIwr_total, 'r');
hold off
grid on

hx = xlabel('$\theta$');
hy = ylabel('$\theta^2 J(\theta)$');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
legend([h1 h2], 'Without Replacement', 'With Replacement', ...
'Location', 'SouthEast');

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/relFIcompare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

figure(2)

h1 = semilogy(theta,1./(theta.^2.*FIwor_total), 'color', dgrey, 'linewidth', linewidth);
hold on
h2 = semilogy(theta,1./(theta.^2.*FIwr_total), 'r');
hold off
grid on

hx = xlabel('$\theta$');
hy = ylabel('CRLB(log $\theta$)');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/relerr_compare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

figure(3)

h1 = semilogy(theta,FIwor_total, 'color', dgrey, 'linewidth', linewidth);
hold on
h2 = semilogy(theta,FIwr_total, 'r');
% h3 = plot(theta,theta.^2.*b./(theta.*(1-theta)),'k');
hold off
grid on
hx = xlabel('$\theta$');
hy = ylabel('$J(\theta)$');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
legend([h1 h2], 'Without Replacement', 'With Replacement', ...
'Location', 'SouthEast');

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/FIcompare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

figure(4)

h1 = semilogy(theta,1./(FIwor_total), 'color', dgrey, 'linewidth', linewidth);
hold on
h2 = semilogy(theta,1./(FIwr_total), 'r');
hold off
grid on

hx = xlabel('$\theta$');
hy = ylabel('CRLB($\theta$)');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/var_compare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

if (save_data)
    outfile = ['EEC_n', num2str(n), '_N', num2str(N), '_b', num2str(b), '_l', num2str(li)];
    save outfile theta FIwor_total FIwr_total
end

% figure(100)
% h1 = semilogy(theta,1./(FIwr_total), 'color', 'r', 'linewidth', linewidth);
% hold on
% h2 = semilogy(theta,1./(FIwr), 'k--', 'linewidth', linewidth);
% hold off
% grid on
% 
% hx = xlabel('$\theta$');
% hy = ylabel('CRLB($\theta$)');
% xlim([lo hi]);
% axes_handle = gca;
% set(axes_handle, 'FontSize', fontsize);
% set(hx, 'FontSize', fontsize);
% set(hy, 'FontSize', fontsize);
