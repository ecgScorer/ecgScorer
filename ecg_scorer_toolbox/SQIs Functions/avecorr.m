% % function m_corr = avecorr(sig,fs,R)
% ------------------------------------------------------------------------
%    avecorr  - Assess the ECG signal Quality by evaluating
%            how the average QRS complex is consistent with each heartbeat
%    Inputs
%       sig : the ECG signal
%       R : position of detected QRS complexes
%       fs : sampling frequency
%
%    Ver. 1.0.0
%  
% ---------------------------------------------------------------------------------------


function m_corr = avecorr(sig,~,R)
  
    %% *********** template matching
   %% limit border

nb_B = length(R);
mRR = min(mean(diff(R)),median(diff(R)));
    med_R = round(mRR);
    beat_on = round(R - (med_R-1)/2); % beginning indice of each the beat
    beat_off = round(R + (med_R-1)/2); % end indice of each beat
%check border limit to avoid negative/non-existant indices
    if beat_on(1)<1
        a=find(beat_on<1,1,'last') +1; % beat indice to start tem0.plate
    else
        a=1;
    end
    if beat_off(nb_B)>length(sig)
        z=find(beat_off>length(sig),1,'first')-1;%nb_B - 1; % beat indice to end template
    else 
        z = nb_B;
    end
%% beats Matrice Generation
%Nbeats = z-a+1;
Beats=zeros(z-a+1,med_R); % number of beat considered by the templates
ii=1;
for i=a:z
    Beats(ii,:) = sig(beat_on(i):beat_off(i));
    ii=ii+1;
end

template = mean(Beats);
    %% 
    % figure, plot(template)
    %correlation coef
    cof=zeros(1,z-a+1);
        for i=1:z-a+1
            cor=corrcoef(template,Beats(i,:));
            cof(i)=cor(1,2);
        end
        m_corr = mean(cof);
        
end