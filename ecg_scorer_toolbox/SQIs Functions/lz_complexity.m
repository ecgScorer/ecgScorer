%% function [complexity, sub_strings]=lz_complexity(sequence)
% ------------------------------------------------------------------------
%    encodingLZC  - Computes the Lempel-Ziv-Welch complexity.
%
%   Input: an array containing a sequence of 0s and 1s.
%   Output: the complexity as the number of words in the dictionary computed 
%           according to the LZW algorithm, and the dictionary of words.
%  
%    Ver. 1.0.0
%  
%    Created:         Fotsing kuetche (23.06.2023)
%                     Noura Alexendre   
%
%                     The University of Ngaoundere
%    mail: ecgScorer@gmail.com
%
%    when using this code cite :
%
%  [1] F. Kuetche, N. Alexendre, N. E. Pascal, and S. Thierry,
%               ‘Simple, efficient, and generalized ECG signal quality assessment method
%               for telemedicine applications’, Inform. Med. Unlocked, 
%               vol. 42, p. 101375, 2023, doi: 10.1016/j.imu.2023.101375.
%
%  [2] Y. Zhang, S. Wei, C. Di Maria, and C. Liu, "Using Lempel-Ziv Complexity to Assess ECG Signal Quality," 
%                J. Med. Biol. Eng., vol. 36, no. 5, pp. 625–634, 2016, doi: 10.1007/s40846-016-0165-5.
%
%============================================================================================================

function [complexity, sub_strings]=lz_complexity(sequence)

if (~isequal(unique(sequence), [0 1]') && ...
    ~isequal(unique(sequence), [0]) && ...
    ~isequal(unique(sequence), [1]))
    error('Input must be an array containing 0s and/or 1s only.');
end
%% converting sequence to string
b='01';
sequence = b(sequence+1);
%% Lempel Ziv complexity

sub_strings ={sequence(1)};
n=length(sequence);
ind = 2;
inc = 0;

 while ind+inc <= n
    
    if ind+inc > n
        break
    end
    
    sub_str = sequence(ind : ind + inc);
    if ismember({sub_str}, sub_strings)
            inc=inc+1;
    else
            sub_strings = [sub_strings(:)', {sub_str}];
            ind =ind+inc+1;
            inc = 0;
    end
end
   complexity = length(sub_strings) ;
   %% normalize
   b = n*log(2)/log(n);
   complexity = complexity/b;
end