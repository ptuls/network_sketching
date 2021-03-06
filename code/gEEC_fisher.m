%% Comparing Fisher information for gEEC scheme for sampling without
%  replacement versus sampling with replacement approximation. Here,
%  immunity was assumed. The Fisher information is computed for any single 
%  level with immunity. 
%
%  Paul Tune Fri 25 Jul 2014
%  paul.tune@adelaide.edu.au
%
%  Last updated: Sat 26 Jul 2014
%  
clear;

% bit flip probability i.e. theta. Range this from 'lo' to 'hi'
lo = 0.001;
hi = 0.6;
lim = 100;
theta = linspace(lo,hi,lim);
save_data = 0;

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
b = 1048;
% Length of sample
li = 32;
% Alphabet size
K = 5;
% Number of non-zero bits
n = 300;
% Packet size in bits (the population in sampling terminology)
N = 1000;

%% How this ties in with previous work:

FIwr_total = zeros(1,lim);
FIwor_total = zeros(1,lim);
comb = zeros(1,b+1);

%% Compute bit sampling probabilities
%  sampling with replacement 
p = n/N;
L = 2^(b*li);
G = b*li;  % number of bits sampled from N
B = 2^b;

% sampling with replacement probabilities (for comparison)
% p1 = [1-p p];
% for i=2:G
%     p1 = kron(p1,[1-p p]);
% end

% sampling without replacement (need to find better way of calculating)
% for j=0:G
%     comb(j+1) = nchoosek(G,j);
% end
% h = hygepdf(0:G,N,n,G)./comb;
% 
% % case input is (0,0,...,0),(0,0,...,1),...,(1,1,...,1) in that order
% p2 = h(sum(de2bi(0:(L-1)),2)+1);  % probability of all input sequences

% need to find the XOR differences between pairs of sequences (transmitter
% vs receiver sequences parity calculations)
%
% note: doing this in C might be much faster, MATLAB bad at bit operations
% and loops
%

% Old MATLAB based code: much much slower
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

% Faster C-based code with bit twiddling
% tic
% table = generate_table(li,b);
% toc

tic

FIwr_level = zeros(lim,li);
for i=1:lim
    m = 1-theta(i);
    md = -1;
    if (K > 1)
        m = [m theta(i)/2];
        md = [md 0.5];
        if (K > 2)
            m = [m zeros(1,K-3) theta(i)/2];
            md = [md zeros(1,K-3) 0.5];
        end
    end
    
    % Construct probabilities and differential
    f = zeros(1,K);
    fd = f;
    f(1) = 1;

    
    for j=1:li
        fd = cconv(fd,m,K) + cconv(f,md,K);
        f = cconv(f,m,K);
        % need to clean from possible convolution artefacts
        f(f < 0) = 0;
        FIwr_level(i,j) = sum(fd(f > 0).^2./f(f > 0));
    end
    
    FIwr_total(i) = b*FIwr_level(i,li);
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

% h1 = semilogy(theta,theta.^2.*FIwor_total, 'color', dgrey, 'linewidth', linewidth);
% hold on
semilogy(theta,theta.^2.*FIwr_level(:,1)', 'r');
hold on
for j=2:li
    semilogy(theta,theta.^2.*FIwr_level(:,j)', 'r');
end
hold off
grid on

hx = xlabel('$\theta$');
hy = ylabel('$\theta^2 J(\theta)$');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
% legend([h1 h2], 'Without Replacement', 'With Replacement', ...
% 'Location', 'NorthWest');

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/relFIcompare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

figure(2)

% h1 = semilogy(theta,1./(theta.^2.*FIwor_total), 'color', dgrey, 'linewidth', linewidth);
% hold on
h2 = semilogy(theta,1./(theta.^2.*FIwr_total), 'r');
% hold off
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

% h1 = semilogy(theta,FIwor_total, 'color', dgrey, 'linewidth', linewidth);
% hold on
h2 = semilogy(theta,FIwr_total, 'r');
% h3 = plot(theta,theta.^2.*b./(theta.*(1-theta)),'k');
% hold off
grid on
hx = xlabel('$\theta$');
hy = ylabel('$J(\theta)$');
xlim([lo hi]);
axes_handle = gca;
set(axes_handle, 'FontSize', fontsize);
set(hx, 'FontSize', fontsize);
set(hy, 'FontSize', fontsize);
% legend([h1 h2], 'Without Replacement', 'With Replacement', ...
% 'Location', 'NorthWest');

if (print_figures)
    set(gcf, 'PaperPosition', position);
    filename = ['Plots/FIcompare_n' num2str(n) '_N' num2str(N) '.eps'];
    print('-depsc', filename);
    fprintf('printing %s\n', filename);
end

figure(4)

% h1 = semilogy(theta,1./(FIwor_total), 'color', dgrey, 'linewidth', linewidth);
% hold on
h2 = semilogy(theta,1./(FIwr_total), 'r');
% hold off
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
