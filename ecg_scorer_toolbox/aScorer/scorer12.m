% ------------------------------------------------------------------------
%    scorer12  - Assess ECG quality using kuetche et al method [1]
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
% function [label, comment, signalNum, aScore, iScore] = scorer12(x,fs)
% Assess the quality for ECG signal Dn
%  
% Inputs:      
%       x: Single or multichannel ECG signal. The channel must be a column vector.
% 
%       fs: the Sampling frequency in Hz
%       
% Outputs:
%       label: The quality class of the ECG. 3 possible values :
%                   0   = Good Quality Signal
%                   0.5 = intermediate Quality signal 
%                   1   = Bad Quality Signal
%       Comment : String describibg the problem of the signal
%               or the number of passage to compute the aScore.
%       signalNum : the order of signal in the provided database 
%       aScore : adjusted ECG quality Score between 0 and 1, 
%       iScore : initial ECG Quality Score between 0 and 1.
%               when iScore/aScore = -1 the signal was detected to be a
%               pure noise or saturated/flatline signal
% Example Usage:
%       [label, comment, signalNum, aScore, iScore] = scorer12(signal,250);
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [label, comment, signalNum, aScore, iScore] = scorer12(x,fs)
Fs = fs;
numberOfSig = size(x,2);
aScore = - ones(numberOfSig,1);
iScore = - ones(numberOfSig,1);
label = rand(numberOfSig,1);
comment = strings(numberOfSig,1);
signalNum = (1:numberOfSig)';
for i=1:numberOfSig
    sig=x(:,i);
%% 1* Flat line detection and saturation detection
flaty = isFlatline(sig,fs);
if flaty ==1 %flat line for 1.5sec
    % existence d'une line plate
    disp('*******************************************************')
    disp("| flat line / saturation / electrode problem detected |")
    disp('*******************************************************')
    class = 1;
    label(i) = class;
    comment(i)="flat line / saturation / electrode problem detected";
    continue
end
%% 2* Pure Noise detection
nzc = isPurenoise(sig,fs);
if nzc>=200
    class = 1;
    label(i) = class;
    disp('***********************')
    disp("| pure noise detected |")
    disp('***********************')
    comment(i) = 'pure noise detected';
    continue
end
if nzc>=160 && nzc< 200
    disp('************************')
    disp("| suspected pure noise |")
    disp('************************')
    comment(i) = 'suspected pure noise';
end

%% pan tompkins detector doubled band pass sqi
sig = normalize(sig,'range',[0, 1]);
[qrs2,i2,~]=pan_tompkin2(sig,fs,8,20);% detector 2

% heuristic rules for bad signals
N2=length(i2);
mm=(N2);
tmSQI = max(diff([1,i2,10*fs]));
flag = 0;

if mm<=3 || tmSQI>2.5*fs
    disp("found less than 4 peaks")
    % search the R peaks using the reversed signal
    disp("---- Search the R peaks using the reversed signal ----")
    sig_rev = sig(end:-1:1);
    [~,irev]=pan_tompkin2(sig_rev,fs,8,20);
    df=([1,irev,10*fs]);
    [maxi, imax] = max(diff(df)); % the max of R-R indices
    if length(irev)<=3 %|| max(diff([1,irev,10*fs]))>2.5*fs
        % 1) the minimal HR = 24 BPM (bradycardia)
        disp('** the number of beats is truly less than 4 **')
        class = 1;
        label(i) = class;
        comment(i) = 'HR found less than 24 BPM';
        continue
    else
        %recover original indices
        i2 = length(sig)+1 - irev;
        i2 = i2(end:-1:1);
        qrs2 = sig(i2);
        flag = 1;
        disp("** mised peaks recovered **")
        %plot(sig), hold on, plot(i2, sig(i2),"*r"), hold off
    end 

    if maxi>=3*fs
        % search missing R peaks stating from the point after the max
        [~,irev2] = pan_tompkin2(sig_rev(df(imax)+10 : end),fs,8,20);
        irev2 = irev2 + df(imax)+10 -1;
        irev = [irev(1:imax-1), irev2]; % add the recovered R peaks
        if max(diff([1,irev,10*fs]))>2.5*fs % if the time is still great
            class = 1;
            label(i) = class;
            comment(i) = "time between two beat surpassed 2.5 sec";
            disp('time between two beat surpassed 2.5 sec')
            continue
        else
            %recover original indices
            i2 = length(sig)+1 - irev;
            i2 = i2(end:-1:1);
            qrs2 = sig(i2);
            flag = 1;
            disp("found r-peaks in the 2.5 sec interval")
            %plot(sig), hold on, plot(i2, sig(i2),"*r"), hold off
        end

    end
end
% disp('computing cons')
%RRi = diff(i2./fs)
%cSQI = std(RRi)/mean(RRi)
% %L24 = length(RRi(RRi>(60/24)));
% %L220 = length(RRi(RRi<(60/220)));
QRS_amp = max(diff(qrs2));


if  QRS_amp>0.8 % ((L24 + L220)>20 && (L24 + L220)<30) ||
    disp('*****************')
    disp('| abrupt Change |')
    disp('*****************')
    class = 1;
    label(i) = class;
    comment(i) = "abrupt Change detected";
    continue
end

%% Beats' average correlation algorithm
L2=0.35;
L3=0.375;
[iscore, ascore]=simSQI(sig,i2,flag,Fs);
iScore(i) = iscore;
aScore(i) = ascore;
%disp('____stage 2___')

if iscore <L2
    class = 1;
    label(i) = class;
    disp('********** Bad Signal *************')
    disp(['the iscore ', num2str(iscore), ' is < ', num2str(L2)])
    continue
elseif iscore>=0.9
    comment(i) = "one pass";
    class = 0;
elseif iscore>=L2 && iscore<0.9 && ascore>=0.9
    comment(i) = "two pass";
    class = 0;
elseif iscore>=L2 && iscore<0.9 && ascore<0.9 && (ascore>=L3 )
    class = 0.5;
else
    class = 1;
end
if class==1
    disp('********** Bad Signal *************')
    disp(['the iscore ', num2str(iscore), ' the aSore is ', num2str(ascore)])
    disp('')
else
    disp('********** Good Signal *************')
    disp(['the iscore ', num2str(iscore), ' the aSore is ', num2str(ascore)])
    disp('')
end
label(i)=class;
end
