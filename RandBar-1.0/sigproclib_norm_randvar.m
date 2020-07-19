
% -----------------------------------------------------------------
%  sigproclib_norm_randvar.m
%
%  This function normalizes a given random variable.
%
%  input:
%  X     - random variable
%
%  output:
%  Xnorm - normalized random variable
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function Xnorm = sigproclib_norm_randvar(X)

    % check number of arguments
    if nargin < 1
        error('Too few inputs.')
    elseif nargin > 1
        error('Too many inputs.')
    end
    
    
    % normalize
    Xnorm = (X - mean(X))/std(X);

return
% -----------------------------------------------------------------
