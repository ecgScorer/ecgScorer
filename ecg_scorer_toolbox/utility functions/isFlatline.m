% function condition = isFlatline(sig,Fs)
% ------------------------------------------------------------------------
%    isFlatline  - Assess if the ECG signal contains flatline or saturation
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


function condition = isFlatline(sig,Fs)
if sum(sig)==0
condition=1; % all zeros, no signal
return
end
 sig = fix(normalize(sig,'range',[0,50]));
%figure,plot(sig)

condition1 = flat_calc(sig,Fs,0);

sigo = sig(sig>(1.5*mean(sig)));
condition11 = flat_calc(sigo,Fs,1);

% sigo = sig(sig<(mean(sig)));
% condition111 = flat_calc(sigo,Fs,1);
condition = or(condition1,(condition11));


function inter = flat_calc(sig,Fs,c)

xx = (diff(sig));
flat_i = (find(xx==0 | xx==-1 | xx==1));
flato = find(diff(flat_i,2)==0);
k=1;
matrix = [];

while k < length(flato)-1
    
    counter =  0;
    while flato(k+1)==flato(k)+1
        counter = counter+1;
        k=k+1;
        if k==length(flato)
        break
        end
    end
    matrix = [matrix, counter];
    k=k+1;
    if k==length(flato)-2
        break
    end
    
end
if c==1
flaty = sum(matrix);
inter = flaty>round(2.5*Fs);
else
flaty = max(matrix);
inter = flaty>round(1.5*Fs); 
end