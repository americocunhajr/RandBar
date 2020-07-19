
% ----------------------------------------------------------------- 
%  randbar_fixed_mass_spring__phys_param.m
%
%  This function defines the simulation parameters for
%  a random 1D fixed-mass-spring bar.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 20, 2012
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function [t0,t1,rho,E,d,A,L,c,k,kNL,m,...
          alpha1,alpha2,sigma,nlflag,...
          Nx,dx,xmesh,dt,time,Ndt,...
          Nmodes,freq_nat,phi,grad_phi,lambda] = ...
          randbar_fixed_mass_spring__phys_param(case_name)

      
% initial time (s)
t0 = 0.0;

% final time (s)
t1 = 20.0e-3;

% mass density (kg/m^3)
rho = 7900.0;

% elastic modulus (Pa)
E = 203.0e9;

% cross section diameter (m)
d = 50.0e-3;

% cross section area (m^2)
A = pi*d^2/4;

% bar length (m)
L = 1.0;

% damping constant (N.s/m)
c = 5.0e3;

% linear stiffness constant (N/m)
k = 650.0;

% nonlinear stiffness constant (N/m^3)
kNL = 650.0e13;
%kNL = 0.0;

% attached mass (kg)
m = 0.1*rho*A*L;

if ( strcmp(case_name,'case1') )
    m = 0.0*rho*A*L;
elseif ( strcmp(case_name,'case1a') )
    m = 0.1*rho*A*L;
elseif ( strcmp(case_name,'case1b') )
    m = 1.0*rho*A*L;
elseif ( strcmp(case_name,'case1c') )
    m = 10.0*rho*A*L;
elseif ( strcmp(case_name,'case1d') )
    m = 50.0*rho*A*L;
end

% initial displacement amplitude 1 (m)
alpha1 = 0.1e-3;

% initial displacement amplitude 2
alpha2 = 0.5e-3;

% external force amplitude (N)
sigma = 5.0e3;
%sigma = 0.0;

% linear (0) / nonlinear (1) flag
nlflag = 1;
%nlflag = 0;



% number of mesh points
Nx = 5*12;

% spatial resolution (m)
dx = L/(Nx-1);

% spatial mesh vector (m)
xmesh = (0.0:dx:L)';



% number of modes
if ( strcmp(case_name,'modes_conv') )
    Nmodes = 30;
else
    Nmodes = 10;
end

% bar modes and natural frequencies
[freq_nat,phi,grad_phi,lambda,~] = ...
randbar_fixed_mass_spring__modes(Nmodes,xmesh,rho,E,A,L,k,m);

% maximum frequency of the band (Hz)
freq_max = max(freq_nat);

% Nyquist frequency (Hz)
freq_nyquist = 2.0*freq_max;

% sampling frequency (Hz)
freq_samp = 10.0*freq_max;

% sample points spacing in the frequency domain (Hz)
dfreq = 50.0;

% time step or sampling period (s)
dt = 1.0/freq_samp;

% initial time (s)
t0 = 0.0;

% final time (s)
t1 = 1.0/dfreq;

% refined temporal mesh vector (s)
time = (t0:dt:t1)';

% number of time steps
Ndt = length(time);

return
% -----------------------------------------------------------------
