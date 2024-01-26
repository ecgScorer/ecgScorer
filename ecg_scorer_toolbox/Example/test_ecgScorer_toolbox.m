%%
% 
%   Created by :         Fotsing kuetche (23.06.2023)
%
%%
% 
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
%  
%    cite as  F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry,
%               ‘Simple, efficient, and generalized ECG signal quality assessment method
%               for telemedicine applications’, Inform. Med. Unlocked, 
%               vol. 42, p. 101375, 2023, doi: 10.1016/j.imu.2023.101375.
%
%% 
% This file show the main functionalities of the software
%% Loading the database
% the 'signal_test' database contains 10 ECG recordings. Those data have been 
% extracted from the Computer in cardiology Challenge 2011 challenge dataset A 
% [21].  It comport 5 signals of goods quality, and 5 of bad quality (2 Gaussian 
% noise, and 3 electrodes problem).

% load the signal
fs = 250; % the sampling frequency of the database is 250 Hz
load('signals_test.mat')

%% 
% _ECGs should be vectors or matrices with leads columnwise arranged (temporal 
% dimension in lines)_
%% Assess the Quality of the Signals
% 1 to 5 are good signals

% ploting the signal
plot(signals(:,3))
xlabel('samples')
ylabel("Amplitude")
% assess the quality
[rep1, com1, sigNum] = scorer12(signals(:,1:5), 250);
% 6 and 7 are gaussian noise

% ploting the signal
plot(signals(:,6))
xlabel('samples')
ylabel("Amplitude")
% assess the quality
[rep2, com2, sigNum2] = scorer12(signals(:,6:7), 250);
% 8 to 10 are saturation

% ploting the signal
plot(signals(:,10))
xlabel('samples')
ylabel("Amplitude")

% assess the quality
rep3 = scorer12(signals(:,8:10), 250);
%% Extract Signal Quality Indices
% Statistics-based SQIs
%%
% 
%  % function output = statSQI(Dn,stat)
%  % Compute Statistics-based Signal Quality Indices specified by stat for ECG signal Dn
%  % 
%  % Inputs:      
%  %   Dn: Single or multichannel ECG signal. The channel must be a column vector.
%  % 
%  %   stat : a string specifying which statistics to output. Possible values are :
%  %             'min' for the minimum, 'mean' for the mean, 'max' for the maximum, 'std' for standard deviation, 
%  %             'var' for variance, 'max-min' the difference between max and min, 
%  %             'max-mean' the difference between max and mean, 
%  %             'MAC' the mean of absolute change, 'LSAM' the longest strike after mean,
%  %             'ASC' the absolute sum of change, 'ST' the sum of the timeseries, 'AE' the absolute energy,
%  %             'nzc' the number of zero crossing, 'skn' the skewness, 'kur' the kurtosis,
%  %             'hosSQI' the hosSQI
%  %             'des', Descriptive statistics (min, mean, max, std, and var)
%  %             'HOS' Higher Order Statistics (skn, kurt, and hosSQI), 
%  %             'all' output all the underlined SQI
%  % Outputs:
%  %  output: The SQIs specified by 'stat'. For multiple SQI and multiple ECG
%  %                     SQIs are row based, while signal is column based.
%
%% 
% 
% 
% Example of execution for all the ECG in the database

statSQI(signals,"mean") % compute the mean off all the ECGs
statSQI(signals,'HOS') % compute the Higher Order Statistics (skn, kurt, and hosSQI)
% Frequency Domain-based SQIs
%%
% 
%  % function output = frequencySQI(Dn,alg, fs)
%  % Compute frequency domain based Signal Quality Indices specified by 'alg', for ECG signal Dn
%  %  
%  % Inputs:      
%  %       Dn: Single or multichannel ECG signal. The channel must be a column vector.
%  % 
%  %       alg : a string specifying which frequency domain-based methods to calculate. Possible values are :
%  % 
%  %       'psdl' for Power Spectral Density (PSD) in low frequencies 0-1Hz
%  %       'psdh' for PSD in high freq 10-100Hz,
%  %       'psdn' for the PSD of the overall signal
%  %       'bassqi' relative power in the baseline
%  %       'iorsqi'  in-band to out-of-band spectral power ratio
%  %       'psqi' power of QRS over the rest of the signal 
%  %       'Lpsqi' power in low freq 0-3Hz over the overall power
%  %       'mpsqi' midle freq 5-35Hz over the overall power
%  %       'hpsqi' power in high freq 40-end Hz over the overall power
%  %       'all' output all the underlined SQIs    
%  %       
%  % Outputs:
%  %       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%  %       SQIs are row based, while signal is column based.
%

% compute the power of QRS over the rest of the signal for the first signal of the database
frequencySQI(signals(:,1),'psqi',fs) 
% compute all the frequency domain based SQIs for the 3rd signal of the database
frequencySQI(signals(:,3),'all',fs) 
% Non linear SQIs
% Function
%%
% 
%  % function output = nonLinearSQI(Dn,alg, varargin)
%  % Compute non-linear Signal Quality Indices specified by 'alg', for ECG signal Dn
%  %  
%  % Inputs:      
%  %       Dn: Single or multichannel ECG signal. The channel must be a column vector.
%  % 
%  %       alg : a string specifying which non-linear methods to output. Possible values are :
%  % 
%  %       'app' for approximate entropy, 'samp' for sample entropy,
%  %       'fuz' for fuzzy entropy, 'fuzM' for fuzzy measure entropy,
%  %       'lzc' for lempel-iv complexity, 'elz' for encoding lzc,
%  %       'h' for hurst exponent 
%  %       'all' output all the underlined SQIs
%  %       
%  %       varargin (optional) : can be one of the variables 'm', 'tau',or 'a'
%  %       followed by the desired value.
%  %       m      --- embedding dimension. default m = 2.
%  %       tau    --- time delay (for sample entropy only). default tau = 1.
%  %       a      --- multiplying coeficient (between 0.15 and 0.25) to compute r, 
%  %             the threshold value to determine similarity r = a*std(Dn). default a = 0.2
%  %       
%  % Outputs:
%  %       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%  %       SQIs are row based, while signal is column based.
%  %  
%  % Example Usage:
%  %       out = nonLinearSQI(Dn,'samp')
%  %       out = nonLinearSQI(Dn,'samp','m',2,'tau',3)
%
% Examples

% compute the sample entropy for 7th signal
%  with default settings (m=2, and tau = 1, a = 0.2)
nonLinearSQI(signals(:,7),'samp')
% The previous statement equivalent (m=2, and tau = 1, a = 0.2)
nonLinearSQI(signals(:,7),'samp','m',2,'tau',1,'a',0.2)
% Now we compute the sample entropy for 7th signal
%  with settings (m=2, and tau = 1, a = 0.18)
nonLinearSQI(signals(:,7),'samp','m',2,'tau',1,'a',0.18)
% QRS-detector based SQIs
% Function
%%
% 
%  % function output = qrsDetectorSQI(Dn,alg, fs)
%  % Compute qrs detector-based Signal Quality Indices specified by 'alg', for ECG signal Dn
%  %  
%  % Inputs:      
%  %       Dn: Single or multichannel ECG signal. The channel must be a column vector.
%  % 
%  %       alg : a string specifying which QRS detector-based methods to calculate. Possible values are :
%  % 
%  %       'tSQI'   average correlation coefficient between beats
%  %       'corSQI' average correlation coefficient between the detected beats and the average beats.
%  %       'iSQI'   the ratio of the 15th percentile value to the 85th percentile value of the RR intervals
%  %       'eSQI'   The relative energy of the QRS complex
%  %       'pcaSQI' the sum of the eigenvalues associated with the five principal components to the sum of all eigenvalues
%  %       'all'    output all the underlined SQIs    
%  %       
%  % Outputs:
%  %       output: The SQI(s) specified by 'alg'. For multiple SQI and multiple ECG
%  %       SQIs are row based, while signal is column based.
%  %  
%  % Example Usage:
%  %       out = qrsDetectorSQI(Dn,'corSQI', 250)
%  %       out = qrsDetectorSQI(Dn, 'tSQI', 1000)
%
% Examples

% compute the corSQI (orphanidou method) for signals 1 and 2
qrsDetectorSQI(signals(:,1:2),'corSQI',fs)
%% 
% For instances involving Gaussian noises and saturated signals, SQIs based 
% on QRS detectors may not be computed, as no QRS complexes are detected. In such 
% cases, the program returns -1 for all SQIs.

% compute all qrs detectors based SQIs for 7th to 10th signal
qrsDetectorSQI(signals(:,7:10),'all',fs)