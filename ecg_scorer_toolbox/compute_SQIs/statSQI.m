% ------------------------------------------------------------------------
%    statSQI  - Compute Statistics-based Signal Quality Indices
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%    
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
% ------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function output = statSQI(Dn,stat)
% Compute Statistics-based Signal Quality Indices specified by stat for ECG signal Dn
%  
% Inputs:      
%       Dn: Single or multichannel ECG signal. The channel must be a column vector.
%
%       stat : a string specifying which statistics to output. Possible values are :
%
%       'min' for the minimum, 'mean' for the mean, 'max' for the maximum, 'std' for standard deviation, 
%       'var' for variance, 'max-min' the difference between max and min, 
%       'max-mean' the difference between max and mean, 
%       'MAC' the mean of absolute change, 'LSAM' the longest strike after mean,
%       'ASC' the absolute sum of change, 'ST' the sum of the timeseries, 'AE' the absolute energy,
%       'nzc' the number of zero crossing, 'skn' the skewness, 'kur' the kurtosis,
%       'hosSQI' the hosSQI
%       'des', Descriptive statistics (min, mean, max, std, and var)
%       'HOS' Higher Order Statistics (skn, kurt, and hosSQI), 
%       'all' output all the underlined SQI
% Outputs:
%       output: The SQIs specified by 'stat'. For multiple SQI and multiple ECG
%       SQIs are row based, while signal is column based.
%  
% Example Usage:
%       out = statSQI(Dn,'all')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = statSQI(Dn,stat)

% Evaluate input argument 'stat'
if sum(strcmp(stat, {'min', 'mean', 'max', 'std', 'var', 'max-min', 'max-mean',...
            'MAC', 'LSAM', 'ASC', 'ST', 'AE', 'nzc', 'skn', ...
            'kur', 'hosSQI', 'des', 'HOS', 'all'}))==0
        error(char([stat, ' is not a recognized argugument', ' you must enter a valid argument']));
end
        
        
% compute basic statistics : min, max, std, and var
mini = min(Dn);
moy = mean(Dn);
maxi = max(Dn);
stde=std(Dn);
varia =var(Dn);
% Difference between max and min/mean
difMax_min = maxi - mini; 
difMax_moy = maxi - moy;
% mean of absolute change MAC
V = diff(Dn); MAC = mean(abs(V));

% absolute sum of change
Abs_sum_chan = sum(abs(V));

% sum over time series
ST = sum(Dn);

% Absolute Energy
AE = sum(Dn.^2);

%Zero crossing rate and Longest Strike After Mean
zci = @(v) find(diff(sign(v)));
nzc=zeros(1,size(Dn,2));

for k = 1:size(Dn,2)

    nzc(1,k)=length(zci(Dn(:,k)));

    %% longest strike above mean
    Dnn=Dn(:,k)- moy(k);
    ind = find((Dnn)>=0);
    Vind = diff(ind);
    vv = find(Vind~=1);
    N=length(Vind);
    if isempty(vv)
        LS= N;
    elseif length(vv)>=2
        for ii=1:length(vv)-1
           if ii==1
               mat(ii)=length(Vind(1:vv(ii)));
           else
               mat(ii)=length(Vind(vv(ii):vv(ii+1)));
           end
        end
        LS = max(mat);
    else
        mat(1)=length(Vind(1:vv(1)));
        mat(2)=length(Vind(vv(1):end));
        LS=max(mat);
    end
    LSAM(1,k)=LS;
end
%%
% Skewness and kurtosis, and higher order stat
skn = skewness(Dn);
kurt = kurtosis(Dn);
hosSQI = abs(skn).*kurt/5;


%convert to line
mini = mini'; moy = moy'; maxi = maxi'; stde = stde'; varia = varia';
 difMax_min = difMax_min'; difMax_moy = difMax_moy';
 MAC = MAC'; LSAM = LSAM'; Abs_sum_chan = Abs_sum_chan'; ST = ST';
 AE =AE'; nzc = nzc'; skn = skn'; kurt = kurt'; hosSQI = hosSQI';
 
% output
    if strcmp(stat,'all')
        output=table(mini, moy, maxi, stde, varia, difMax_min, difMax_moy,...
            MAC, LSAM, Abs_sum_chan, ST, AE, nzc, skn, kurt, hosSQI);
    elseif strcmp(stat,'min')==1
        output = table(mini);
    elseif strcmp(stat,'mean')==1
        output = table(moy);
    elseif strcmp(stat,'max')==1
        output = table(maxi);
    elseif strcmp(stat,'std')==1
        output = table(stde);
    elseif strcmp(stat,'var')==1
        output = table(varia);
    elseif strcmp(stat,'max-min')==1
        output = table(difMax_min);
    elseif strcmp(stat,'max-mean')==1
        output = table(difMax_mean);
    elseif strcmp(stat,'MAC')
        output = table(MAC);
    elseif strcmp(stat,'LSAM')
        output = table(LSAM)';
    elseif strcmp(stat,'ASC')
        output = (Abs_sum_chan)';
    elseif strcmp(stat,'ST')==1
        output = table(ST);
    elseif strcmp(stat,'AE')
        output = table(AE);
    elseif strcmp(stat,'nzc')
        output = table(nzc);
    elseif strcmp(stat,'skn')
        output = table(skn);
    elseif strcmp(stat,'kur')
        output = table(kurt);
    elseif strcmp(stat,'hosSQI')==1
        output = table(hosSQI);
    elseif strcmp(stat,'des') % Descriptives statistics
        output = table(mini, moy, maxi, stde, varia);
    elseif strcmp(stat,'HOS') %higher order statistics
        output = table(skn, kurt, hosSQI);
    else
        disp('you must specify a correct argument')
    end

end
