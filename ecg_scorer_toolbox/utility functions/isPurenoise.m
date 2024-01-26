% function nzc = isPurenoise(v,fs)
% ------------------------------------------------------------------------
%    isPurenoise  - Compute the Number of Zeros Crossing (nzc) for the given
%    signal v, having a sampling frequency fs
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (20.06.2023)
%    
%                     The University of Ngaoundere
%    mail: fotsing.fk@gmail.com
%    cite as [1] F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry,
%               ‘Simple, efficient, and generalized ECG signal quality assessment method
%               for telemedicine applications’, Inform. Med. Unlocked, 
%               vol. 42, p. 101375, 2023, doi: 10.1016/j.imu.2023.101375.
% ---------------------------------------------------------------------------------------

function nzc = isPurenoise(v,fs)
[b,a] = butter(3,0.5./(fs/2),'high');
v = filtfilt(b,a,v);
v = normalize(v,"range",[-1,1]);
zci = find(diff(sign(v)));
nzc = length(zci);
end
