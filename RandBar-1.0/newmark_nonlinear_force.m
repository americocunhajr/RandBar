
% -----------------------------------------------------------------
%  newmark_nonlinear_force.m
%
%  This function computes nonlinear external force matrix.
%
%  Input:
%  Udisp      - Ndofs x 1 displacement vector
%  phi        - Nx x Ndofs modal matrix
%  Nx         - number of mesh points
%  phys_param - physical parameters vector
%
%  Output:
%  FNL        - Ndofs x Ndt nonlinear external force
%  dFNL_dU    - Ndofs x Ndofs nonlinear external force derivative matrix
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 13, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [FNL,dFNL_dU] = newmark_nonlinear_force(Udisp,phi,Nx,phys_param)
    
    % check number of arguments
    if nargin < 4
        error(' Too few inputs.')
    elseif nargin > 4
        error(' Too many inputs.')
    end
    
     % physical parameters
%     rho    = phys_param(1);        % mass density
%     E      = phys_param(2);        % elastic modulus
%     A      = phys_param(3);        % cross section area
%     L      = phys_param(4);        % bar length
%     k      = phys_param(5);        % linear stiffness constant
     kNL    = phys_param(6);        % nonlinear stiffness constant
%     alpha1 = phys_param(7);        % initial displacement amplitude 1
%     alpha2 = phys_param(8);        % initial displacement amplitude 2
%     sigma  = phys_param(9);        % external force amplitude
%     c      = phys_param(10);       % damping constant
%     m      = phys_param(11);       % attached particle mass
    
    
    % compute nonlinear external force
    Udisp_L = phi(Nx,:)*Udisp;
    phi_L   = phi(Nx,:)';
    FNL     = -kNL*phi_L*(Udisp_L.^3);
    
    % compute nonlinear external force derivative
    phi_ij  = phi(Nx,:)'*phi(Nx,:);
    dFNL_dU = -3*kNL*phi_ij*(Udisp_L.^2);
    
return
% -----------------------------------------------------------------