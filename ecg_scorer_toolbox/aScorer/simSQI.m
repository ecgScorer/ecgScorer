% ----------------------------------------------------------------------------------------
%    simSQI  - compute the SQIs iScore and aScore 
%               proposed by kuetche et al method [1]
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%    
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
%    cite as [1] F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry,
%               ‘Simple, efficient, and generalized ECG signal quality assessment method
%               for telemedicine applications’, Inform. Med. Unlocked, 
%               vol. 42, p. 101375, 2023, doi: 10.1016/j.imu.2023.101375.
% ---------------------------------------------------------------------------------------
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% [iScore, aScore]= simSQI(sig,i2,flag,fs)
% 
% Inputs:      
%       x: Single-channel ECG signal. The channel must be a column vector.
%       i2: detected R-peaks indices matrix
%       flag: 1 or 0, the R-peaks indices matrix is provided or not
%       fs: the Sampling frequency in Hz
%       
% Outputs:
%       aScore : adjusted ECG quality Score between 0 and 1, 
%       iScore : initial ECG Quality Score between 0 and 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [iScore, aScore]= simSQI(sig,i2,flag,fs)
 %fs = 200;
%% load and display the signal

if flag == 1
    %disp("flag on")
    R = i2;
    %figure,plot(sig), hold on, plot(R, sig(R),'*r'), hold off
else

    [b,a] = butter(3,[1 40]./(fs/2),'bandpass');
    sig = filtfilt(b,a,sig);
    [a,R,b]=pan_tompkin2(sig,fs,8,20);
    if length(R)<4
    [~,R,~]=pan_tompkin2(sig,fs,5,15);% detector 2
    end
end
%figure,plot(sig), hold on, plot(R, sig(R),'*r'), hold off
[iScore, aScore]=runSimSQI(sig,R);

end    
    function [iScore, aScore]=runSimSQI(sig,R)
%% limit border
%Pnoise =[1 0.8 0.6 0.4 0.2 0.025 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%Pgroup =[1 0.6 0.4 0.2 0.1 0.05 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
Pgroup = @(x) 1.7*exp(-0.52.*x) -0.005;
Pnoise = @(x) 1.6*exp(-0.38.*x) -0.01;

% 
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
Nbeats = z-a+1;
Beats=zeros(z-a+1,med_R); % number of beat considered by the templates
ii=1;
for i=a:z
    Beats(ii,:) = sig(beat_on(i):beat_off(i));
    ii=ii+1;
end

%% correlation matrices generation

corMat=corrcoef(Beats');
corMat_mean = mean(corMat);
iScore = mean(corMat_mean);

%% Similarity check algorithm
thr1 = 0.65;
thr2 = 0.75;
        %find the maximum corr
[spB ,idx_maxi] = max(corMat_mean);
        % start research in correlation matrix
idx_worst = find(corMat(idx_maxi,:)<thr1); 
idx_good  = find(corMat(idx_maxi,:)>=thr1 & corMat(idx_maxi,:)~=1);
        % verify if they are simillar cc=corMat(idx_worst,idx_worst)
L=length(idx_worst); 
%clear score_wg;
score_wg = [];
% score_n=[];
if L==1
%     disp('possible noise')
    noise_count=1;
%     cc=(corMat_mean(idx_worst(1)));
%         score_wg = cc*Pnoise(noise_count);
    group=0;
else
        count=1;
        noise_count=0;
        group=0;
    while (count<=L)||(isempty(idx_worst))==1
        if (isempty(idx_worst))==1
            break
        end

        idx_res = corMat(idx_worst(1),idx_worst)>thr2;
        idx_res = idx_worst(idx_res) ;
        if length(idx_res)==1
            noise_count = noise_count+1;
            cc=(corMat_mean(idx_res(1)));
            score_wg(count)= cc*Pnoise(noise_count);
            idx_worst = idx_worst(2:end);
            count=count+1;
            continue
        end
        % comparing indices
        [comp, loc] = ismember(idx_worst,idx_res);
        b_res = (find(comp==1));
        group=group+1;
        sw =mean(mean(corMat(idx_worst(b_res),idx_worst(b_res))));
        score_wg(count) = sw*Pgroup(group);
        idx_worst =idx_worst(comp==0); %update the list
        count=count+1;

    end
end
%% *** final score

Ngood=length(idx_good) +1;
% vScore = (Ngood - noise_count)/Nbeats;

if Ngood==1
    score_gg=spB;
elseif Ngood==2
     score_gg=(spB+corMat(idx_maxi,idx_good))/2;
else
 score_gg = mean(mean(corMat([idx_maxi,idx_good],[idx_maxi,idx_good])));
end

% group;
if ~isempty(score_wg)
    aScore = (mean(score_wg) + score_gg)/2 ;
else
    aScore = score_gg ;
end

end
