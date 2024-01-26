% ------------------------------------------------------------------------
%    qrsDetectorSQI  - Compute QRS detectors-based Signal Quality Indices
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%    
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
% ------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% function output = qrsDetectorSQI(Dn,alg, fs)
% Compute qrs detector-based Signal Quality Indices specified by 'alg', for ECG signal Dn
%  
% Inputs:      
%       Dn: Single or multichannel ECG signal. The channel must be a column vector.
% 
%       alg : a string specifying which qrs-detectors based methods to calculate. Possible values are :
% 
%       'tSQI'   average correlation coefficient between beats
%       'corSQI' average correlation coefficient between the detected beats and the average beats.
%       'iSQI'   the ratio of the 15th percentile value to the 85th percentile value of the RR intervals
%       'eSQI'   The relative energy of the QRS complex
%       'pcaSQI' the sum of the eigenvalues associated with the five principal components 
%                       to the sum of all eigenvalues
%       'all'    output all the underlined SQIs    
%       
% Outputs:
%       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%       SQIs are row based, while signal is column based.
%  
% Example Usage:
%       out = qrsDetectorSQI(Dn,'corSQI', 250)
%       out = qrsDetectorSQI(Dn, 'tSQI', 1000)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
function feat = qrsDetectorSQI(X,alg,Fs)

ind = strcmp(alg, {'tSQI', 'corSQI', 'iSQI', 'eSQI', 'pcaSQI', 'all'}) ;
ind = find(ind,1);

if sum(strcmp(alg, {'tSQI', 'corSQI', 'iSQI', 'eSQI', 'pcaSQI', 'all'}))==0
        error([alg, ' is not a recognized argument', ' you must enter a valid argument'])
end

[~, N] = size(X);
logic = floor(1*heaviside(-ind+10) + 18*heaviside(ind-10)); %gives 1 (if ind<9), 9 (if ind=9)
feat = zeros(logic, N); % initialize matrix, size = 1*N or 9*N, all = 9 SQIs

for i = 1:N % compute SQI for each channel
     x = X(:,i);

    [~,R_tom,~]=pan_tompkin2(x,Fs,5,15);% detector 2
    N_tom = length(R_tom) ;
    t1 = floor(0.07*Fs);
            t2 =floor( 0.08*Fs);
            lx = length(x);
    
    if N_tom<4
        feat(:,i) = -ones(1,logic);
        continue
    end

    switch alg 
        case 'tSQI'
          %2 tSQI
            feat(:,i) = t_SQI(x,R_tom);
        case 'corSQI'
            %3 average template corSQI
            feat(:,i) = avecorr(x,Fs,R_tom);
        case 'iSQI'
             %4 iSQI
             RR = quantile(R_tom,[0.15,0.85]);
            feat(:,i) = RR(1)/RR(2);
        case 'eSQI'
            %6 eSQI
            
            Er =0;
            for k=1:N_tom
                i1=max(1,R_tom(k)-t1);
                i2=min(lx,R_tom(k)+t2);
            Er = Er + sum((x(i1:i2)).^2);
            end
        
            feat(:,i) = Er/(sum(x.^2));
        case 'pcaSQI'
            %7 pcaSQI
                %time aligne QRS
            t1 = floor(0.1*Fs);
            t2 =t1;
            t_aligned_ECG =[];
            for k=1:N_tom
                i1=(R_tom(k)-t1);
                i2=(R_tom(k)+t2);
                if i1<1
                    continue
                elseif i2>lx
                    break
                end
            t_aligned_ECG = [t_aligned_ECG,x(i1:i2)];
            end
            [~,~,latent] = pca(t_aligned_ECG);
           
            
            ss = size(t_aligned_ECG,2);
            if ss<=5 && ss>3
                aa = 3;
            elseif ss<=3
                aa = 1;
            else
                aa=5;
            end
            
            feat(:,i) = sum(latent(1:aa))/sum(latent);
        case 'all'
            tSQI = t_SQI(x,R_tom);

                %3 average template corSQI
                corSQI = avecorr(x,Fs,R_tom);
                
                %4 iSQI
                RR = quantile(R_tom,[0.15,0.85]);
            %     RR85 = quantile(R_tom,0.85);
                iSQI = RR(1)/RR(2);
            
             
                %6 eSQI
                t1 = floor(0.07*Fs);
                t2 =floor( 0.08*Fs);
                lx = length(x);
                Er =0;
                for k=1:N_tom
                    i1=max(1,R_tom(k)-t1);
                    i2=min(lx,R_tom(k)+t2);
                Er = Er + sum((x(i1:i2)).^2);
                end
            
                eSQI = Er/(sum(x.^2));
            
                %7 pcaSQI
                    %time aligne QRS
                t1 = floor(0.1*Fs);
                t2 =t1;
                t_aligned_ECG =[];
                for k=1:N_tom
                    i1=(R_tom(k)-t1);
                    i2=(R_tom(k)+t2);
                    if i1<1
                        continue
                    elseif i2>lx
                        break
                    end
                t_aligned_ECG = [t_aligned_ECG,x(i1:i2)];
                end
                [~,~,latent] = pca(t_aligned_ECG);
               
                
                ss = size(t_aligned_ECG,2);
                if ss<=5 && ss>3
                    aa = 3;
                elseif ss<=3
                    aa = 1;
                else
                    aa=5;
                end
                
                pcaSQI = sum(latent(1:aa))/sum(latent);
            

            feat(1:5,i) = [tSQI, corSQI, iSQI, eSQI, pcaSQI]';
            name = {'tSQI', 'corSQI', 'iSQI', 'eSQI', 'pcaSQI'};
    end
end
    
if strcmp(alg,'all') == 0
    name = {alg};
end

feat = array2table(feat','VariableNames',name);
