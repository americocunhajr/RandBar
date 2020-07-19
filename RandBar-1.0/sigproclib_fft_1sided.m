
% -----------------------------------------------------------------
%  sigproclib_fft_1sided.m
%
%  This function computes the onesided FFT of a signal x(t) in
%  time domain and returns the corresponding transformned signal
%  X(freq) in frequency domain.
%
%  input:
%  x    - signal in time domain
%  fs   - sampling rate (Hz)
%  nfft - FFT number of points (optional)
%
%  output:
%  X    - transformed signal in frequency domain
%  freq - onesided frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [X,freq] = sigproclib_fft_1sided(x,fs,nfft)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 3
        error('Too many inputs.')
    elseif nargin == 2
        nfft = length(x);
    end
    
    % check arguments
    if fs <= 0.0
        error('fs must positive.')
    end
    
    
    % sampling time
    T = nfft*(1/fs);
    
    
    % frequency sampling points
    k = 0:1:nfft-1;
    
    
	% frequency vector
    freq = k/T;
    
    
	% compute FFT
	X = fft(x,nfft);
    
    
    % normalize FFT data
    X = X/sqrt(nfft);
    
    
    % take only the first half of the spectrum
    cut_off = ceil(nfft/2) + 1;
    X       = X(1:cut_off);
    freq    = freq(1:cut_off);

return
% -----------------------------------------------------------------
