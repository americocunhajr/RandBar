
% -----------------------------------------------------------------
%  randbar_fixed_mass_spring__modes.m
%
%  This function computes natural frequencies and modes of a
%  fixed-mass-spring bar.
%  
%  Input:
%  Nmodes - number of modes desired
%  x      - position vector
%  rho    - mass density (kg/m^3)
%  E      - elastic modulus (Pa)
%  A      - cross section area (m^2)
%  L      - bar length (m)
%  k      - spring stiffness constant (N/m)
%  m      - punctual mass (kg)
%
%  Output:
%  freq_nat - natural frequencies vector (Hz)
%  phi      - modes matrix (Nx x Nmodes)
%  grad_phi - modes derivatives matrix (Nx x Nmodes)
%  lambda   - caracteristic parameters vector
%  ising    - number of singularities found (optional)
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 20, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [freq_nat,phi,grad_phi,lambda,ising] = ...
         randbar_fixed_mass_spring__modes(Nmodes,x,rho,E,A,L,k,m)
    
    % check number of arguments
    if nargin < 8
        error('Too few inputs.')
    elseif nargin > 8
        error('Too many inputs.')
    end
    
    % error tolerance
    tol = 1.0e-7;
    
    % check input arguments
    if ( rho < 0.0 || abs(A - 0.0) < tol)
        error('rho must be positive.')
    end

    if ( A < 0.0 || abs(A - 0.0) < tol)
        error('A must be positive.')
    end

    if ( L < 0.0 || abs(L - 0.0) < tol)
        error('L must be positive.')
    end    
    
    if ( E < 0.0 )
        error('E must be nonnegative.')
    end
    
    if ( k < 0.0 )
        error('k must be nonnegative.')
    end

    if ( m < 0.0 )
        error('m must be nonnegative.')
    end

    % x vector dimension
    Nx  = length(x);
    
    % preallocate memory for caracteristic parameters vector
    lambda = zeros(Nmodes,1);
    
    % preallocate memory for natural frequencies vector
    freq_nat = zeros(Nmodes,1);
    
    % preallocate memory for modes matrix
    phi = zeros(Nx,Nmodes);

    % preallocate memory for modes derivatives matrix
    grad_phi = zeros(Nx,Nmodes);
    
    % root counter
    iroot = 0;
    
    % singularity counter
    ising = 0;
    
    % bar wave speed
    c = sqrt(E/rho);
    
    % caracteristic function
    F = @(lamb) lamb.*cot(lamb) + ((k*L)/(E*A)) - (lamb.^2).*(m/(rho*A*L));
    
    % caracteristic parameters
    dlamb = 1.0e-4*Nmodes;
    lamb1 = 0.0;
    lamb2 = 0.0;
    
    
    while( iroot < Nmodes )
        
        % caracteristic parameter
        lamb1 = lamb1 + dlamb;
        lamb2 = lamb1 + dlamb;
        
        % caracteristic function at lamb_n1 and lamb_n2
        F1 = F(lamb1);
        F2 = F(lamb2);
        
        % checking the necessary condition to be a root
        if F1*F2 < 0.0
            
            % looking for a root or a singularity of F
            [lamb_n,Fn] = fzero(F,[lamb1,lamb2]);
            
            % check if lamb_n found is a root or a singularity
            if ( abs(Fn) < tol )
                
                iroot             = iroot + 1;
                lambda(iroot,1)   = lamb_n;
                freq_nat(iroot,1) = lamb_n*(c/L);
                phi(:,iroot)      = sin(lamb_n*(x/L));
                grad_phi(:,iroot) = (lamb_n/L)*cos(lamb_n*(x/L));
                
            else
                
                ising = ising + 1;
                
            end
            
        end
        
    end
    
    % convert from (rad/s) to Hz
    freq_nat = freq_nat/(2*pi);

return
% -----------------------------------------------------------
