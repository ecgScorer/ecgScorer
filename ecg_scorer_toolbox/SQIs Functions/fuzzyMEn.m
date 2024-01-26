function y = fuzzyMEn(x,m)
% calculation of the fuzzy Measure Entropy
% X = {x1, x2, ..., xn} is a sequence
% m is the embedding dimension. m=2 (default)
% m<N-1

%  m=2
% x=mit_nor(1,:);%rand(1,2500);%[1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2 1 2 2];
if nargin==1
m = 2;
end
[phi_Lm, phi_Fm] = loc_glob_fuz(x,m);
[phi_Lm1, phi_Fm1] = loc_glob_fuz(x,m+1);
%% Fuzzy local/global measure entropy and fuzzyMEn (y)
fL = log(phi_Lm) - log(phi_Lm1);
fg = log(phi_Fm) - log(phi_Fm1);
y = fL + fg;
% entr = fuzzyEn(x, m, 2, 0.2.*std(x))
% FuzEn = FuzzyEn2(x,m,0.2.*std(x),2)
return

function [phi_Lm, phi_Fm] = loc_glob_fuz(x,m)
%% (1) generate Xm segments

N = length(x);
XLm = zeros(m,N-m+1);
XFm = XLm;
for i=1:N-m+1
    val = x(i:i+m-1);
XLm(1:m,i) = val' - mean(val);
XFm(1:m,i) = val' - mean(x);
end
%% (2) distance of local sequence
dLm = zeros(N-m+1,N-m+1);
dFm = dLm;
for i=1:N-m+1
    for j=1:N-m+1
       dLm(i,j) = max(abs(XLm(:,i) - XLm(:,j)));
       dFm(i,j) = max(abs(XFm(:,i) - XFm(:,j)));
    end
end

%% local fuzzy function and global fuzzy function
% rF and rL denote the thresholds of gaussian function
% The nL and nF are their weights of sequence segments similarity
rL = 0.2.*std(x); rF = rL;
nL = 2;
nF = 3;
DLm = exp(-(dLm./rL).^nL);
DFm = exp(-(dFm./rF).^nF);
phi_Lm =0;
phi_Fm = 0;
for k = 1:N-m
    for j=1:N-m
        if j ==k
            continue
        end
        phi_Lm = phi_Lm + DLm(k,j);
        phi_Fm = phi_Fm + DFm(k,j);
    end
end
div = ((N-m)*(N-m-1));
phi_Lm = phi_Lm./div;
phi_Fm = phi_Fm./div;
return