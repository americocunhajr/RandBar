
% -----------------------------------------------------------------
%  randbar_fixed_mass_spring__galerkin.m
%
%  This function computes Galerkin matrices and vectors.
%  
%  Input:
%  Nbasis     - number of basis functions
%  time       - time vector
%  lambda     - lambda parameter vector
%  Fwn        - white noise external force
%  phys_param - physical parameters vector
%
%  Output:
%  M          - mass matrix
%  C          - damping matrix
%  K          - stiffness matrix
%  F          - force vector
%  U0         - initial displacment vector
%  V0         - initial velocity vector
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 20, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [M,C,K,F,U0,V0] = ...
randbar_fixed_mass_spring__galerkin(Nbasis,time,lambda,Fwn,phys_param)
    
    % check number of arguments
    if nargin < 5
        error('Too few inputs.')
    elseif nargin > 5
        error('Too many inputs.')
    end
    
    % check input parameter
    if ( Nbasis <= 0 )
        error('Nbasis must be a positive integer.')
    end
    
    % physical parameters
    rho    = phys_param(1);        % mass density
    E      = phys_param(2);        % elastic modulus
    A      = phys_param(3);        % cross section area
    L      = phys_param(4);        % bar length
    k      = phys_param(5);        % linear stiffness constant
%    kNL    = phys_param(6);        % nonlinear stiffness constant
    alpha1 = phys_param(7);        % initial displacement amplitude 1
    alpha2 = phys_param(8);        % initial displacement amplitude 2
    sigma  = phys_param(9);        % external force amplitude
    c      = phys_param(10);       % damping constant
    m      = phys_param(11);       % attached particle mass
    
    % wave speed
    cL = sqrt(E/rho);
    
    % number of time steps
    Ndt = length(time);
    
	% preallocate memory for mass matrix
	M = zeros(Nbasis,Nbasis);
    
	% preallocate memory for damping matrix
	C = zeros(Nbasis,Nbasis);
    
	% preallocate memory for stiffness matrix
	K = zeros(Nbasis,Nbasis);
    
	% preallocate memory for force vector
	F = zeros(Nbasis,Ndt);
    
    % preallocate memory for initial displacement vector
	U0 = zeros(Nbasis,1);
    
    % preallocate memory for initial velocity vector
	V0 = zeros(Nbasis,1);
    
    
    % compute Galerkin matrices and vectors
    l1    = lambda(1);
    c_l1  = cos(l1);
    s_l1  = sin(l1);
    
    for i=1:1:Nbasis
    
    	li   = lambda(i);
    	s_li = sin(li);
    	c_li = cos(li);
        
        for j=1:1:Nbasis
            
            % mass matrix
            M_ii   = 0.5*rho*A*L*(li - s_li*c_li)/li + m*s_li^2;
            M(i,j) = KronD(i,j)*M_ii;

            % damping matrix
            C_ii   = 0.5*c*L*(li - s_li*c_li)/li;
            C(i,j) = KronD(i,j)*C_ii;
            
            % stiffness matrix
            K_ii   = 0.5*(E/L)*A*li*(c_li*s_li + li) + k*s_li^2;
            K(i,j) = KronD(i,j)*K_ii;
            
        end
        
        % excitation frequency
        omega = l1*(cL/L);
        
        % external force
        F1     = 0.5*sigma*L*Fwn(:)*(l1 - c_l1*s_l1)/l1;
        F(i,:) = F1*KronD(i,1);
%        F(i,:) = 0.0;  sin(omega*time)*
        
        % initial displacement
        U0(i,1) = alpha1*KronD(i,3) + ...
                  2*alpha2*L*(li*c_li - s_li)/(li*c_li*s_li - li^2);
        
        % initial velocity
        V0(i,1) = 0.0;
        
    end
    

    
return
% -----------------------------------------------------------------
