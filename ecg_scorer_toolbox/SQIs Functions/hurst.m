% Hurst exponent
% Author: Conrado Chiarello
% Date: 02.01.2019
% NUEM - Multiphase flow research center
% UTFPR - Technological University of Parana
% modified version 25.08.2023
% by Fotsing Kuetche

function H = hurst(x)
% only accept timeseries signals x = [a1, a2, a3, ..., aN]

x=x(:)'; %force to timeseries. vectorize, then transpose 

% Creates local matrix for calculations
X = tril(repmat(x, length(x), 1));
size(X);
X(1,:) = [];

% Makes zeros equal NaN for calculations
X(X == 0) = NaN;

% Range and Std Dev
meanMat = repmat(mean(X,2,'omitnan'),1,length(x));
R = max(cumsum(X -meanMat,2),[],2) - min(cumsum(X - meanMat,2),[],2);
S = nanstd(X,0,2);

% % Plot for visualization
%plot(log10(2:length(x)), log10((R./S)'), '-')

% Power law fitting
powerLaw = polyfit(log10(2:length(x)), log10((R./S)'), 1);

% Hurst exponent
H = powerLaw(1);

