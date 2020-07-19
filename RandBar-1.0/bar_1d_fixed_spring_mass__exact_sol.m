
% -----------------------------------------------------------------
%  bar_1D_solution.m
%
%  This function computes the exact solution for the
%  fixed-mass-spring bar problem.
%
%  Input:
%  x      - position vector (m)
%  time   - time vector (s)
%  Nmodes - number of modes in the series
%
%  Output:
%  u      - displacement field
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: May 18, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function u = rand_1d_fixed_spring_mass_bar__exact_sol(x,time,Nmodes,phys_param)
    
    % check number of arguments
    if nargin < 4
        error(' Too few inputs.')
    elseif nargin > 4
        error(' Too many inputs.')
    end
    
    % dimension of time vector
    Ndt = length(time);
        
    % preallocate memory for a(t)
    a = zeros(Nmodes,Ndt);

    % preallocate memory for Cn
    Cn = zeros(Nmodes,1);

    % preallocate memory for Sn
    Sn = zeros(Nmodes,1);
    
    % physical parameters
    rho    = phys_param(1);        % mass density
    E      = phys_param(2);        % elastic modulus
    A      = phys_param(3);        % cross section area
    L      = phys_param(4);        % bar length
    k1     = phys_param(5);        % linear stiffness constant
%    k2    = phys_param(6);        % nonlinear stiffness constant
    alpha1 = phys_param(7);        % initial displacement amplitude 1
    alpha2 = phys_param(8);        % initial displacement amplitude 2
%    nu1    = phys_param(9);        % first natural frequency
%    sigma  = phys_param(10);       % external force amplitude
%    gamma  = phys_param(11);       % natural frequency scalling factor
    c      = phys_param(12);       % damping constant
    m      = phys_param(13);       % attached particle mass
    
    % wave speed
    cbar = sqrt(E/rho);
    
    % damping factor
    ksi = c/(E*A);
    
    % compute bar nomral modes parameters
    % -----------------------------------------------------------
    [~,phi,~,lambda] = bar_1D_modes(Nmodes,x,rho,E,A,L,k1,m);
    % -----------------------------------------------------------
    
    
    for i=1:1:Nmodes
        
        % compute lambda
        li   = lambda(i);
        c_li = cos(li);
        s_li = sin(li);
        
        % damped frequency
        wd = sqrt((li/L)^2 - 0.25*(cbar*ksi)^2)*cbar;
        
        % cumpute Cn coeficient
        Cn(i,1) = alpha1*KronD(i,3) + ...
                  2*alpha2*L*(li*c_li - s_li)/(li*c_li*s_li - li^2);
        
        % cumpute Sn coeficient
        Sn(i,1) = 0.5*((ksi/wd)*cbar^2)*Cn(i,1) + 0.0;
                       
        % cumpute time dependent coeficients
        a(i,:)  = exp(-0.5*ksi*cbar^2*time).*(Cn(i,1)*cos(wd*time) + ...
                  Sn(i,1)*sin(wd*time));
    end
    
    % compute u(x,t)
    u = phi*a;
    
return
% -----------------------------------------------------------------
