
% -----------------------------------------------------------------
%  sigproclib_psd_1sided.m
%
%  This function computes the onesided power spectral density (PSD)
%  of a given signal x(t) in time domain.
%
%  input:
%  X    - (1 x Ndt) signal in time domain
%  fs   - sampling rate (Hz)
%  nfft - FFT number of points (optional)
%  win  - (1 x Ndt) window (optional)
%
%  output:
%  PSD_X - onesided power spectral density of X
%  freq  - onesided frequency vector (Hz)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Mar 6, 2014
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [PSD_X,freq] = sigproclib_psd_1sided(X,fs,nfft,win)
    
    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 4
        error('Too many inputs.')
    elseif nargin == 2
    	nfft = length(X);
        win  = ones(nfft);
    elseif nargin == 3
        win = ones(length(X));
    end
    
    % check arguments
    if fs <= 0.0
        error('fs must positive')
    end
    
    % convert X to a row vector (if necessary)
    if find( size(X) == max(size(X)) ) < 2
        X = X';
    end
    
    % convert win to a row vector (if necessary)
    if find( size(win) == max(size(win)) ) < 2
        win = win';
    end
    
    
    % time step (s)
    dt = 1/fs;
    
    
    % compute FFT and frequency vector
    [Xfft,freq] = sigproclib_fft_1sided(X.*win,fs,nfft);
	
    
    % compute PSD
    PSD_X = dt*abs(Xfft).^2;
    
    
    % conservation of total power
    %(zero frequency and the Nyquist frequency do not occur twice.)
    PSD_X(2:end-1) = 2*PSD_X(2:end-1);
	
return
% -----------------------------------------------------------------
