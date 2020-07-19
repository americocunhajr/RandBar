
% -----------------------------------------------------------------
%  norm_H1.m
%
%  This function computes the H1 norm of a diferentiable function
%  u: [xmin,xmax] -> R with finite support.
%
%                   -- xmax         -- xmax
%                  |               |
%  normL2 := sqrt( | u(x)^2 dx +   | u'(x)^2 dx )
%                  |               |
%                -- xmin         -- xmin
%  
%  Input:
%  xmin   -   low limit of integration
%  xmax   - upper limit of integration
%  u      - function u
%  du     - function u derivative
%
%  Output:
%  normL2 - norm L2
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Oct 31, 2011
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function normH1 = norm_H1(xmin,xmax,u,du)
    
    % check number of arguments
    if nargin < 4
        error(' Too few inputs.')
    elseif nargin > 4
        error(' Too many inputs.')
    end
    
    % check input arguments
    if ( xmax < xmin )
        error('xmax must be greather than xmin.')
    end
    
    % number of points in u vector
    N1 = length(u);
    
    % number of points in du vector
    N2 = length(du);
    
    if ( N1 ~= N2 )
        error(' vectos u and du must have the same dimension.')
    end
    
    dx = (xmax-xmin)/(N1-1);
    
    % compute independent variable vector
    x = xmin:dx:xmax;
    
    % compute square od function u
    u2 = u.*u;
    
    % compute square od function du
    du2 = du.*du;
    
    % compute H1 norm
    normH1 = sqrt(trapz(x,u2) + trapz(x,du2));
    
return
% -----------------------------------------------------------------
