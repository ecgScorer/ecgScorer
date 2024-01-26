function  [yf, tf]=spectre(f,fs)
% that return the power spectral density (yf)
% and the frequence interval (tf)
% f input signal. if f is a matrix, spectre treats each column
% fs sampling frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f=f(:)';
l =length(f);
yf=fft(f);
yf_norm =yf.*conj(yf)./l ;
yf_norm=yf_norm(1:floor(l/2)) ;
tf =fs*(0:l/2-1)/l ;
yf = yf_norm ;

end