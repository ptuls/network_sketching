%% Comparing Fisher information for EEC scheme for sampling without
%  replacement versus sampling with replacement approximation. The Fisher
%  information is computed for group 1, i.e. the parity bit equals to the
%  sampled bit, with immunity. Here, the sketch size was set at b = 2 bits. 
%  This can be generalised further for the other levels and b > 2, but 
%  there are insights already from this simple model.
%
%  Paul Tune Tue 15 Jul 2014
%  paul.tune@adelaide.edu.au
%
%  Last updated: Wed 23 Jul 2014
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
% Model: packet contents -> sampled bits -> channel -> parity bits
% (equations)
%
% Compute Fisher information of received parity bits.

% Group size
b = 2;
% Level
li = 6;
% Number of non-zero bits
n = 100;
% Packet size in bits (the population in sampling terminology)
N = 1000;

%% Interesting points (for b = 2):
% n/N = 0,1: this represents the case when the packet sent has all 0 
% or 1 bits. It is equivalent to a kind of repetition code. In this case,
% sampling with replacement and without replacement yields the same amount 
% of Fisher information. Coincidentally, the Fisher information is equal to
% b/(theta*(1-theta)), which is the bound in "Towards Optimal
% Error-Estimating Codes through the Lens of Fisher Information Analysis".
% Moreover, when n/N is small (or correspondingly, n/N close to 1) there is
% little difference between sampling with and without replacement.
%
% n/N = 0.5: this is the worst case input for both sampling with and
% without replacement. Here, sampling with replacement gives a CRLB of
% infinity, while sampling without replacement gives a large CRLB,
% especially at theta = 0.25. The reason for the former is that with
% sampling with replacement, because the proportion of 0 and 1 bits are
% equally likely, estimator is unable to tell if the received bit has been
% flipped or not. Sampling with replacement does better because of the
% extra information of knowing how many 1 bits are left, but even so, at
% theta = 0.25, the Fisher information almost becomes degenerate as the
% input probabilities are almost uniform.
%

%% How this ties in with previous work:
% The Fisher information derived in "Towards Optimal Error-Estimating Codes 
% through the Lens of Fisher Information Analysis" doesn't take into
% account the packet's content when deriving the bounds. There the bound
% derived is equal to the best case scenario when the packet's bits are all
% 0s or 1s. 
%
% What's interesting here is the interplay between the packet contents and
% the bound. When the entropy of the packet content is highest, i.e. when
% n/N = 0.5, for sampling with replacement, the Fisher information becomes
% zero. The lesson here is to eliminate sources of variability, so the best
% strategy is to send a packet with all 0s or all 1s so what you get are
% just the bit flips from the channel. 

FIwr_total = zeros(1,lim);
FIwor_total = zeros(1,lim);
comb = zeros(1,b+1);

%% Compute bit sampling probabilities
%  sampling with replacement 
p = n/N;
L = 2^(b*li);
G = b*li;  % number of bits sampled from N
B = 2^b;

% sampling without replacement (need to find better way of calculating)
for j=0:G
    comb(j+1) = nchoosek(G,j);
end
h = hygepdf(0:G,N,n,G)./comb;

% case input is (0,0,...,0),(0,0,...,1),...,(1,1,...,1) in that order
pin2 = h(sum(de2bi(0:(L-1)),2)+1);  % probability of all input sequences


% need to find the XOR differences between pairs of sequences (transmitter
% vs receiver sequences parity calculations)
%
% note: doing this in C might be much faster, Matlab bad at bit operations
% and loops
%
% tic
% table = zeros(L,L);
% for i=1:L
%     for j=1:L
%         v = bi2de(xor(de2bi(i-1,L),de2bi(j-1,L)));
%         table(i,j) = bi2de(selfxor(v,G,b));
%     end
% end
% toc

tic
table = generate_table(li,b);
toc

tic
Td1 = [-1 1; 1 -1];
for i=1:lim
    %% Fisher information: sampling with replacement 
    t = [1-theta(i) theta(i)];

    %% Fisher information: sampling without replacement
    T1 = toeplitz(t',t);
    T2 = T1;
    for k=1:(G-1)
        T2 = kron(T2,T1);
    end   
    T2 = diag(pin2)*T2;
    
    % calculating differential w.r.t. theta
    Td2 = zeros(L,L);
    for j=1:G
        if (j==1)
            P = Td1;
        else
            P = T1;
        end
        for k=2:G
            if (j==k)
                P = kron(P,Td1);
            else
                P = kron(P,T1);
            end
        end
        Td2 = Td2 + P;
    end
    Td2 = diag(pin2)*Td2;
    
    f2 = zeros(1,B);
    fd2 = zeros(1,B);
    for j=0:(B-1)
        f2(j+1) = sum(T2(table == j));
        fd2(j+1) = sum(Td2(table == j));
    end

%     for j=1:L
%         for k=1:L
%             v = bi2de(xor(de2bi(j-1,L),de2bi(k-1,L)));
%             id = bi2de(selfxor(v,G,b));
%             f2(id+1) = f2(id+1) + T2(j,k);
%             fd2(id+1) = fd2(id+1) + Td2(j,k);
%         end
%     end
    
    FIwor_total(i) = sum(fd2.^2./f2);
    % since i.i.d., so multiply with subsketch size 
    FIwr_total(i) = b*(4*li^2*(1-2*theta(i))^(2*li-2))/(1-(1-2*theta(i))^(2*li));
end
toc

print_figures = 0;
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
'Location', 'NorthWest');

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
'Location', 'NorthWest');

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
