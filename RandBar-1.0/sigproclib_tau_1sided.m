
% -----------------------------------------------------------------
%  sigproclib_tau_1sided.m
%
%  This function computes time lag vector.
%
%  input:
%  fs   - sampling rate (Hz)
%  nfft - FFT number of points
%
%  output:
%  tau - time lag vector (s)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function tau = sigproclib_tau_1sided(fs,nfft)

    % check number of arguments
    if nargin < 2
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    end
    
    
    % onesided FFT number of points
    nfft = nfft/2;
    
    % sampling points vector
    if mod(nfft,2) == 0
        % if nfft is even
        k = -nfft/2:1:nfft/2-1;
    else
        % if nfft is odd
        k = -(nfft-1)/2:1:(nfft-1)/2;
    end
    
    % new time step
    dt = 1.0/fs;
    
    % time lag vector
    tau = k*dt;
    
return
% -----------------------------------------------------------------
