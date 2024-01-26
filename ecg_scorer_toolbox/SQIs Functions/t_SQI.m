% % function [score]= t_SQI(sig,R)
% ------------------------------------------------------------------------
%    t_SQI  - Assess the ECG signal Quality by evaluating the morphological
%            coherence of the signal over time.
%    Inputs
%       sig : the ECG signal
%       R : position of detected QRS complexes
%
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (01.01.2023)
%                     Noura Alexendre
%
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
% ---------------------------------------------------------------------------------------

function [score]= t_SQI(sig,R)

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

%% correlation matrices generation

corMat=corrcoef(Beats');
corMat_mean = mean(corMat);
score = mean(corMat_mean);
