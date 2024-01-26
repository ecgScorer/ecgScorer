% ------------------------------------------------------------------------
%    nonLinearSQI  - Compute non-linear Signal Quality Indices
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%    
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
% ------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function output = nonLinearSQI(Dn,alg, varargin)
% Compute non-linear Signal Quality Indices specified by 'alg', for ECG signal Dn
%  
% Inputs:      
%       Dn: Single or multichannel ECG signal. The channel must be a column vector.
% 
%       alg : a string specifying which non-linear methods to output. Possible values are :
% 
%       'app' for approximate entropy, 'samp' for sample entropy,
%       'fuz' for fuzzy entropy, 'fuzM' for fuzzy measure entropy,
%       'lzc' for lempel-iv complexity, 'elz' for encoding lzc,
%       'h' for hurst exponent 
%       'all' output all the underlined SQIs
%       
%       varargin (optional) : can be one of the variables 'm', 'tau',or 'a'
%       followed by the desired value.
%       m      --- embedding dimension. default m = 2.
%       tau    --- time delay (for sample entropy only). default tau = 1.
%       a      --- multiplying coeficient (between 0.15 and 0.25) to compute r, 
%             the threshold value to determine similarity r = a*std(Dn). default a = 0.2
%       
% Outputs:
%       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%       SQIs are row based, while signal is column based.
%  
% Example Usage:
%       out = nonLinearSQI(Dn,'samp')
%       out = nonLinearSQI(Dn,'samp','m',2,'tau',3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function feat =nonLinearSQI(X,alg,varargin)% non linear SQI extraction

ind = strcmp(alg, {'app', 'samp', 'fuz', 'fuzM', 'lzc', 'elz', 'h', 'all'});
ind = find(ind,1);
if sum(strcmp(alg, {'app', 'samp', 'fuz', 'fuzM', 'lzc', 'elz', 'h', 'all'}))==0
        error([alg, ' is not a recognized argument', ' you must enter a valid argument'])
end

% initialize default values
m=2; tau = 1; a = 0.2;
% inputs = {'m', 'tau', 'a'};
%overide default vlues if any entries

for n=1:2:nargin-2
    if(~isempty(varargin{n}))
        switch varargin{n}
            case 'm'
                m = varargin{n+1};
            case 'tau'
                tau = varargin{n+1};
            case 'a'
                a = varargin{n+1};
        end
   
    end
end

%get dimension of the signal
[~, N] = size(X);
logic = floor(1*heaviside(-ind+8) + 14*heaviside(ind-8)); %gives 1 (if ind<8), 8 (if ind=8)
feat = zeros(logic, N); % initialize matrix, size = 1*N or 7*N
%tic
parfor i = 1:N % compute SQI for each channel

    % Default values
    x = X(:,i);
    r = a*std(x);
    n = 2;
    % coarse graining process for lzc and elz
    xg = coarseGraining(x);


    switch ind
        case 1
        % Approximate entropy
        feat(:,i) = approximateEntropy(x);
        case 2
        % sample entropy
        feat(:,i) = sampEn(x,m,tau,r);
        case 3
        % fuzzy Entropy
        feat(:,i) = FuzzyEn2(x,m,r,n);
        case 4
        % fuzzy Measure entropy
        feat(:,i) = fuzzyMEn(x,m);
        case 5
        % Lempel-ziv complexity
        feat(:,i) = lz_complexity(xg);
        case 6
        % Encoding Lempel-ziv
        feat(:,i) = encodingLZC(x);
        case 7
        % hurst Exponent
        feat(:,i) = hurst(x);

        case 8
        app = approximateEntropy(x);
        samp = sampEn(x,m,tau,r);
        fuz = FuzzyEn2(x,m,r,n);
        fuzM = fuzzyMEn(x,m);
        lz_c = lz_complexity(xg);
        elz = encodingLZC(x);
        h = hurst(x);
        feat(:,i) = [app, samp, fuz, fuzM, lz_c, elz, h]';
       
    end

end
%toc

if strcmp(alg,'all') == 0
    name = {alg};
else
     name = {'app', 'samp', 'fuz', 'fuzM', 'lzc', 'elz', 'h'};
end

feat = array2table(feat','VariableNames',name);

end