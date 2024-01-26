% function complexity = encodingLZC(x)
% ---------------------------------------------------------------------------------------
%    encodingLZC  - Compute the encoding lempel-ziv complexity
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%                     Noura Alexendre
%
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
%    
% whe using this code cite 
%
%  [1] F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry,
%               ‘Simple, efficient, and generalized ECG signal quality assessment method
%               for telemedicine applications’, Inform. Med. Unlocked, 
%               vol. 42, p. 101375, 2023, doi: 10.1016/j.imu.2023.101375.
%
%  [2] Y. Zhang, S. Wei, H. Liu, L. Zhao, and C. Liu,
%                “A novel encoding Lempel–Ziv complexity algorithm for quantifying
%                 the irregularity of physiological time series,” 
%                 Comput. Methods Programs Biomed., vol. 133, pp. 7–15,
%                 Sep. 2016, doi: 10.1016/j.cmpb.2016.05.010.
%=========================================================================================
function complexity = encodingLZC(x)
% x is matrix where each collumm represents a signal (vector) [----]';

n = length(x); %length of the signal
N = 3*n;        % length of the binarise signal
%% The first binary digit b1(i)
b1 = coarseGraining(x);
%% The second binary digit b2(i)
dx = diff(x);
thr=zeros(1,size(x,2));
bin_dx = coarseGraining(dx,thr);
b2 = zeros(size(x));
b2(2:end,:) = bin_dx;
%% the third binary digit b3(i)
dm = mean(dx);
flag = coarseGraining(dx,dm);
b3 = zeros(size(x));
b3(2:end,:) = not(xor(bin_dx,flag));

%% recontruct the final output
y = zeros(N,size(x,2));
%indices of b1, b2, and b3 in the final signal
i1 = 1:3:N;
i2 = 2:3:N;
i3 = 3:3:N;
y(i1,:) = b1(:,:);
y(i2,:) = b2(:,:);
y(i3,:) = b3(:,:);
%% compute the complexity
 complexity = lz_complexity(y);
return
