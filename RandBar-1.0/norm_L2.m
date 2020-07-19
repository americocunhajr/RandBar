
% -----------------------------------------------------------------
%  norm_L2.m
%
%  This function computes the L2 norm of a function
%  u: [xmin,xmax] -> R with finite support.
%
%                   -- xmax
%                  |
%  normL2 := sqrt( | u(x)^2 dx )
%                  |
%                -- xmin
%  
%  Input:
%  xmin   - low   limit of integration
%  xmax   - upper limit of integration
%  u      - function u to be integrated
%
%  Output:
%  normL2 - L2 norm of u
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 1, 2011
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function normL2 = norm_L2(xmin,xmax,u)
    
    % check number of arguments
    if nargin < 3
        error(' Too few inputs.')
    elseif nargin > 3
        error(' Too many inputs.')
    end
    
    % check input arguments
    if ( xmax < xmin )
        error('xmax must be greather than xmin.')
    end
    
    % number of points in u vector
    N = length(u);
    
    % compute x vector step
    dx = (xmax-xmin)/(N-1);
    
    % compute independent variable vector
    x = xmin:dx:xmax;
    
    % compute square function
    u2 = u.*u;
    
    % compute L2 norm
    normL2 = sqrt(trapz(x,u2));
    
return
% -----------------------------------------------------------------
