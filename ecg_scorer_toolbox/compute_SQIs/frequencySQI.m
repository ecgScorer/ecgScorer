% ------------------------------------------------------------------------
%    frequencySQI  - Compute frequency domain based Signal Quality Indices
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%                     Noura Alexendre
%    
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
% ------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function output = frequencySQI(Dn,alg, fs)
% Compute frequency domain based Signal Quality Indices specified by 'alg', for ECG signal Dn
%  
% Inputs:      
%       Dn: Single or multichannel ECG signal. The channel must be a column vector.
% 
%       alg : a string specifying which frequency domain-based methods to calculate. Possible values are :
% 
%       'psdl' for Power Spectral Density (PSD) in low frequencies 0-1Hz
%       'psdh' for PSD in high freq 10-100Hz,
%       'psdn' for the PSD of the overall signal
%       'bassqi' relative power in the baseline
%       'iorsqi'  in-band to out-of-band spectral power ratio
%       'psqi' power of QRS over the rest of the signal 
%       'Lpsqi' power in low freq 0-3Hz over the overall power
%       'mpsqi' midle freq 5-35Hz over the overall power
%       'hpsqi' power in high freq 40-end Hz over the overall power
%       'all' output all the underlined SQIs    
%       
% Outputs:
%       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%       SQIs are row based, while signal is column based.
%  
% Example Usage:
%       out = frequencySQI(Dn,'psdn',1000)
%       out = frequencySQI(Dn,'iorsqi',360)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function feat = frequencySQI(X,alg,Fs)

ind = strcmp(alg, {'psdl', 'psdh', 'psdn', 'bassqi', 'iorsqi', 'psqi',...
        'lpsqi', 'mpsqi', 'hpsqi', 'all'}) ;
ind = find(ind,1);

if sum(strcmp(alg, {'psdl', 'psdh', 'psdn', 'bassqi', 'iorsqi', 'psqi',...
        'lpsqi', 'mpsqi', 'hpsqi', 'all'}))==0
        error([alg, ' is not a recognized argument', ' you must enter a valid argument'])
end

[~, N] = size(X);
logic = floor(1*heaviside(-ind+10) + 18*heaviside(ind-10)); %gives 1 (if ind<9), 9 (if ind=9)
feat = zeros(logic, N); % initialize matrix, size = 1*N or 9*N, all = 9 SQIs

for i = 1:N % compute SQI for each channel
    x = X(:,i);
    [fecg, ~] = spectre(x,Fs);
    A = Fs/size(x,1);
    %desired freuencies in Hz
    f0=0;
    f1=1;
    f3=3;
    f5=5;
    f10=10; 
    f15=15; 
    f35=35;
    f40=40;
    %indices of frequencies
    i0 =round(1+ f0/A); i1=round(1+ f1/A); i40=round(1+ f40/A);
    i5=round(1+ f5/A); i3 = round(1+ f3/A);i10=round(1+ f10/A); 
    i15=round(1+ f15/A); i35=round(1+ f35/A);

    switch alg
        case 'psdl'
            % PSD_l 0-1Hz
            feat(:,i) = sum(fecg(i0:i1));
        case 'psdh'
            % PSD_h, 10-100
            feat(:,i) = sum(fecg(i10:end));
        case 'psdn'
            % PSD_n
            feat(:,i) = sum(fecg);
        case 'bassqi'
            % relative power in the baseline
             feat(:,i)=sum(fecg(i1:i40))/sum(fecg(1:i40));
        case 'iorsqi'
            % in-band to out-of-band spectral power ratio
            feat(:,i) = sum(fecg(i5:i40))/(sum(fecg) - sum(fecg(i5:i40)));
        case 'psqi'
            % pSQI power of QRS over the rest
            feat(:,i) = sum(fecg(i5:i15))/sum(fecg(i5:i40));
        case 'lpsqi'
            % LpSQI
            feat(:,i) = sum(fecg(i0:i3))/sum(fecg);
        case 'mpsqi'
            % MpSQI
           feat(:,i) = sum(fecg(i5:i35))/sum(fecg);
        case 'hpsqi'
            % HpSQI
           feat(:,i) = sum(fecg(i40:end))/sum(fecg);
           
        case 'all'
        % PSD_l 0-1Hz
        PSDl = sum(fecg(i0:i1));
        % PSD_h, 10-100
        PSDh = sum(fecg(i10:end));
        % PSD_n
        PSDn = sum(fecg);
        % % PSD_h/n
        % rat_h_n = PSDh./PSDn;
        % % PSD_l/n
        % rat_l_n = PSDl./PSDn;

        % relative power in the baseline
            basSQI=sum(fecg(i1:i40))/sum(fecg(1:i40));

        % in-band to out-of-band spectral power ratio
            iorSQI = sum(fecg(i5:i40))/(sum(fecg) - sum(fecg(i5:i40)));

        % pSQI power of QRS over the rest
            pSQI = sum(fecg(i5:i15))/sum(fecg(i5:i40));
        % LpSQI
            LpSQI = sum(fecg(i0:i3))/sum(fecg);
        % MpSQI
            MpSQI = sum(fecg(i5:i35))/sum(fecg);
        % HpSQI
            HpSQI = sum(fecg(i40:end))/sum(fecg);

        feat(1:9,i) = [PSDl, PSDh, PSDn, basSQI, iorSQI, pSQI, LpSQI, MpSQI, HpSQI]';
        name = {'psdl', 'psdh', 'psdn', 'bassqi', 'iorsqi', 'psqi',...
        'lpsqi', 'mpsqi', 'hpsqi'};
    end
end

if strcmp(alg,'all') == 0
    name = {alg};
end
feat = array2table(feat','VariableNames',name);
  
 