
% -----------------------------------------------------------------
%  sigproclib_freq_1sided.m
%
%  This function computes onesided frequency vector.
%
%  input:
%  fs   - sampling rate (Hz)
%  nfft - FFT number of points
%
%  output:
%  freq - onesided frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function freq = sigproclib_freq_1sided(fs,nfft)

    % check number of arguments
    if nargin < 1
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    
    % sampling time
    T = nfft/fs;
    
    % frequency sampling points
    k = 0:1:nfft-1;
    
    % frequency vector
    freq = k/T;
    
    % take only the first half of the spectrum
    cut_off = ceil(nfft/2) + 1;
    freq    = freq(1:cut_off);

return
% -----------------------------------------------------------------
