function entr = sampEn(ts, m, tau, r)
%SAMPEN Sample entropy calculation
%       
%       Inputs:
%               ts     --- time-series
%               m      --- embedding dimension
%               tau    --- time delay
%               r      --- threshold value to determine similarity
% 
%       Ref.:
%       [1] Richman JS, Randall Moorman J. "Physiological time-series
%           analysis using approximate entropy and sample entropy", Am J
%           Physiol Heart Circ Physiol 278: H2039-H2049. 2000.
%       [2] Li P, Liu CY, Li K, Zheng D, Liu C, Hou Y. "Assessing the
%           complexity of short-term heartbeat interval series by
%           distribution entropy", Med Biol Eng Comput 53(1): 77-87. 2015.
% 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%                      (C) Peng Li 2015 -
% If you use the code, please make sure that you cite reference [1]-[2]
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%
% $Author:      Peng Li
%               School of Control Science and Engineering,
%                   Shandong University, PR China
%                   Email: pli@sdu.edu.cn
% $Date:        Mar 23, 2015
% $Modif.:      
% 

% parse inputs
% narginchk(4, 4);

% normalization
ts    = zscore(ts(:));

% reconstruction
N    = length(ts);
indm = hankel(1:N-m*tau, N-m*tau:N-tau);    % indexing elements for dim-m
inda = hankel(1:N-m*tau, N-m*tau:N);        % for dim-m+1
indm = indm(:, 1:tau:end);
inda = inda(:, 1:tau:end);
ym   = ts(indm);
ya   = ts(inda);

% using pdist for saving time, but a large RAM is required
cheb = pdist(ym, 'chebychev'); % inf-norm
cm   = sum(cheb <= r)*2 / (size(ym, 1)*(size(ym, 1)-1));

cheb = pdist(ya, 'chebychev');
ca   = sum(cheb <= r)*2 / (size(ya, 1)*(size(ya, 1)-1));

entr = -log(sum(ca) / sum(cm));