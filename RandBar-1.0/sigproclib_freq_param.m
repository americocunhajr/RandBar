
% -----------------------------------------------------------------
%  sigproclib_freq_param.m
%
%  This function computes parameters necessary to signal analysis
%  in frequency domain given time parameters.
%
%  input:
%  dt   - time step
%  Ndt  - number of time steps
%  Tss  - steady state instant estimation
%
%  output:
%  fs   - sampling frequency (Hz)
%  fcut - cutoff frequency (Hz)
%  nfft - FFT number of points
%  freq - onesided frequency vector (Hz)
%  win  - window vector
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 13, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [fs,fcut,nfft,freq,win] = sigproclib_freq_param(dt,Ndt,Tss)

% sampling frequency (Hz)
fs = 1.0/dt;

% cutoff frequency (Hz)
fcut = 0.5*fs;

% zero Pad factor
zeroPad = nextpow2(10*Ndt);

% number of FFT points
nfft = 2^(zeroPad);

% frequency range vector (Hz)
freq = sigproclib_freq_1sided(fs,nfft);

% window vector
win = window(@tukeywin,Ndt+1-Tss,0.1);

return
% -----------------------------------------------------------------
