
% -----------------------------------------------------------------
%  sigproclib_time_step.m
%
%  This function computes the appropriate time step for a
%  problem with natural frequencies in wn vector.
%  
%  Input:
%  wn      - natural frequencies vector (rad/s)
%  epsilon - time mesh adjust parameter (optional)
%
%  Output:
%  dt      - time step (s)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 7, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function dt = sigproclib_time_step(wn,epsilon)
    
    % check number of arguments
    if nargin < 1
        error('Too few inputs.')
    elseif nargin > 2
        error('Too many inputs.')
    elseif nargin == 1
        epsilon = 0.5;
    elseif ( epsilon >= 1.0 )
        error('dt_param must be less than 1.')
    end
    
    
    % time vector (s)
    tn  = 2*pi./wn;
    
    % time step
    dt = abs(epsilon)*min(tn);
    
return
% -----------------------------------------------------------
