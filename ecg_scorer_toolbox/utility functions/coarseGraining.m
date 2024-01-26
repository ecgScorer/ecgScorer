function y = coarseGraining(x,thr)
% x is a matrix where each columm is a signal
% y is the coarse-grained output
% thr is the imposed threshold, length of thr should be
% equal to number of signal
if nargin==1
thr  = mean(x);    
end
y = ones(size(x));
for i=1:size(x,2)
    ind = x(:,i)<thr(1,i);
    y(ind,i) = 0;
end