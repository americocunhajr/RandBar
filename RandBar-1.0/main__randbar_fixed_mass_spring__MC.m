
% ----------------------------------------------------------------- 
%  main__randbar_fixed_mass_spring__MC.m
%
%                                 |--> u(x,t,w)
%                                  ___
% /|______________________________|   |      k     |\
% /|                              |   |---/\/\/\---|\
% /| rho, E, A                    | m |            |\
% /|______________________________|   |---/\/\/\---|\
% /|                              |___|      kNL   |\
%
%  |----> x
%
%  rho A u_tt =  [E A u_x]_x + f        for 0 < x < L
%           u =  0                      for x = 0
%     E A u_x = -k.u -kNL.u^3 - m.u_tt  for x = L
%           u =  u0                     for t = 0
%         u_t =  v0                     for t = 0
%
%  This script is the main file for a program that simulates
%  the stochastic dynamics of a fixed-mass-spring bar using 
%  Monte Carlo method.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Nov 20, 2012
% -----------------------------------------------------------------


clc
clear all
close all



% program header
% -----------------------------------------------------------
disp('                                                    ')
disp(' ---------------------------------------------------')
disp(' Random Fixed-Spring-Mass Bar                       ')
disp('                                                    ')
disp('                                 |--> u(x,t,w)      ')
disp('                                  ___               ')
disp(' /|______________________________|   |      k     |\')
disp(' /|                              |   |---/\/\/\---|\')
disp(' /| rho, E, A                    | m |            |\')
disp(' /|______________________________|   |---/\/\/\---|\')
disp(' /|                              |___|      kNL   |\')
disp('                                                    ')
disp('  |----> x                                          ')
disp('                                                    ')
disp('  rho A u_tt =  [E A u_x]_x + f        for 0 < x < L')
disp('           u =  0                      for x = 0    ')
disp('     E A u_x = -k.u -kNL.u^3 - m.u_tt  for x = L    ')
disp('           u =  u0                     for t = 0    ')
disp('         u_t =  v0                     for t = 0    ')
disp('                                                    ')
disp('                                                    ')
disp(' by                                                 ')
disp(' Americo Barbosa da Cunha Junior                    ')
disp(' americo.cunhajr@gmail.com                          ')
disp(' ---------------------------------------------------')
disp('                                                    ')
% -----------------------------------------------------------




% name of the case that will be simulated
% -----------------------------------------------------------
%case_name = 'case1';
%case_name = 'case1a';
%case_name = 'case1b';
%case_name = 'case1c';
case_name = 'case1d';


disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------


% simulation parameters
% -----------------------------------------------------------

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
c = 5.0e3; %(mass study)
%c = 5.0; %(NL spring study)

% linear stiffness constant (N/m)
k = 650.0;

% nonlinear stiffness constant (N/m^3)
kNL = 650.0e13;
%kNL = 0.0;

% attached mass (kg)
m = 0.1*rho*A*L;

if ( strcmp(case_name,'case1') )
    m = 0.0*rho*A*L;
    %kNL = 0.0;
elseif ( strcmp(case_name,'case1a') )
    m = 0.1*rho*A*L;
    %kNL = 650.0e12;
elseif ( strcmp(case_name,'case1b') )
    m = 1.0*rho*A*L;
    %kNL = 650.0e13;
elseif ( strcmp(case_name,'case1c') )
    m = 10.0*rho*A*L;
    %kNL = 650.0e14;
elseif ( strcmp(case_name,'case1d') )
    m = 50.0*rho*A*L;
    %kNL = 650.0e15;
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

% PSD estimator upper bound
%eps_v = 0.04;

% number of independent copies of the signal
Ncopies = 1;

% number of time steps in a signal copy
Nsamp = 462;

% minimun frequency of the band (Hz)
freq_min = 0.0;

% maximum frequency of the band (Hz)
freq_max = max(freq_nat);

% sampling frequency (Hz)
freq_samp = 2.0*freq_max;

% frequency resolution or sampling frequency step (Hz)
dfreq = freq_samp/Nsamp;

% time step or sampling period (s)
dt = 1.0/freq_samp;

% period of acquisition for a signal copy (s)
T_samp = Nsamp*dt;

% initial time of analysis (s)
t0 = 0.0;

% final time of analysis (s)
t1 = (Ncopies+1)*T_samp;

% refined temporal mesh vector (s)
time = (t0:dt:t1-dt)';

% number of time steps
Ndt = length(time);

% zero Pad factor
zeroPad = nextpow2(10*Ndt);

% number of FFT points
nfft = 2^(zeroPad);

% steady state instant estimation
if Ncopies == 0
    Tss = 1;
else
    Tss = Nsamp+1;
end

% positions of interest
X1    = floor((3/3)*Nx);
X2    = floor((2/3)*Nx);

% instants of interest
T1    = floor((2/2)*Ndt);
T2    = floor((1/2)*Ndt);

% frequency range vector (Hz)
freq = sigproclib_freq_1sided(freq_samp,nfft);

% window vector
win = window(@tukeywin,Ndt+1-Tss,0.1);
% -----------------------------------------------------------



% pre processing
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Pre Processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% RNG seed
rng_stream = RandStream('mt19937ar','Seed',30081984);
RandStream.setDefaultStream(rng_stream);

% number of samples
Ns = 4^6;
%Ns = 4^4;

if Ns > 1
    
    % elastic modulus mean
    mu_E = E;
    
    % elastic modulus dispersion
    delta_E = 0.1;
    
    % random elastic modulus (Pa)
    E = gamrnd(1.0/delta_E^2,mu_E*delta_E^2,[Ns,1]);
    E(Ns,1) = mu_E;
    
    % nonlinear stiffness mean
    %mu_kNL = kNL;
    
    % nonlinear stiffness dispersion
    %delta_kNL = 0.2;
    
    % random nonlinear stiffness (N/m^3)
    %kNL = gamrnd(1.0/delta_kNL^2,kNL*delta_kNL^2,[Ns,1]);
    %kNL(Ns,1) = mu_kNL;
end

% percentile
ptile = 99;

% cutoff frequency (Hz)
fcut = 0.5*freq_samp;

% white noise spectral density amplitude
S0 = 1.0;

% white gaussian noise standard deviation
sigma_WN = sqrt(4.0*pi*fcut*S0);

% normalized gaussian white noise
WN       = normrnd(0.0,1.0,[Ns Ndt]);
WN(Ns,:) = ones(1,Ndt);


% preallocate memory for u(X,t,w)
U_x1 = zeros(Ns,Ndt);
%U_x2 = zeros(Ns,Ndt);

% preallocate memory for ut(X,t,w)
Ut_x1 = zeros(Ns,Ndt);
%Ut_x2 = zeros(Ns,Ndt);

% preallocate memory for ux(X,t,w)
% Ux_x1 = zeros(Ns,Ndt);
% Ux_x2 = zeros(Ns,Ndt);

% preallocate memory for u(x,T,w)
U_t1 = zeros(Ns,Nx);
% U_t2 = zeros(Ns,Nx);

% preallocate memory for ut(x,T,w)
Ut_t1 = zeros(Ns,Nx);
% Ut_t2 = zeros(Ns,Nx);

% preallocate memory for ux(x,T,w)
Ux_t1 = zeros(Ns,Nx);
% Ux_t2 = zeros(Ns,Nx);

% preallocate memory for FFT of u(X,t,w)
 U_x1_fft = zeros(1,nfft/2+1);
% U_x2_fft = zeros(1,nfft/2+1);

% preallocate memory for FFT of ut(X,t,w)
Ut_x1_fft = zeros(1,nfft/2+1);
%Ut_x2_fft = zeros(1,nfft/2+1);

% preallocate memory for spectral density of u(X,t,w)
PSD_U_x1 = zeros(1,nfft/2+1);
%PSD_U_x2 = zeros(1,nfft/2+1);

% preallocate memory for spectral density of ut(X,t,w)
PSD_Ut_x1 = zeros(1,nfft/2+1);
%PSD_Ut_x2 = zeros(1,nfft/2+1);

% iteration counter vector
iter_vec = zeros(Ns,1);

% vector to store the L2 norm of the 2-norm of u(t)
normL2_Q = zeros(1,Ns);

toc
% -----------------------------------------------------------



% Monte Carlo simulation
% -----------------------------------------------------------
tic

disp(' ');
disp(' --- Monte Carlo Simulation --- ');
disp(' ');

% parallel simulation parameters
%Nprocessors = 8;

%if Ns > Nprocessors^2
%    parflag = 1;
%else
    parflag = 0;
%end


% start MATLAB parallel pool
if parflag == 1
    disp(' ');
    matlabpool('open',Nprocessors)
    disp(' ');
end


parfor_progress(Ns);
for imc=1:Ns
    
    % physical parameters vector
    %phys_param = [rho E(imc) A L k kNL(imc) alpha1 alpha2 sigma c m];
    %phys_param = [rho E A L k kNL(imc) alpha1 alpha2 sigma c m];
    phys_param = [rho E(imc) A L k kNL alpha1 alpha2 sigma c m];

    % Galerkin vectors and matrices
    [M,C,K,F,U0,V0] = ...
     randbar_fixed_mass_spring__galerkin(Nmodes,time,lambda,...
                                           WN(imc,:),phys_param);
    
    % numerical integration of Galerkin ODE system
    [Uj,Utj,Uttj,iter_vec(imc)] = ...
     newmark_ode_solver(M,C,K,F,U0,V0,...
                        dt,Ndt,Nx,Nmodes,phi,nlflag,phys_param);

    % compute U, Ut, and Ux from coeficients and modes
    U_x1(imc,:) = phi(X1,:)*Uj;
%    U_x2(imc,:) = phi(X2,:)*Uj;
     
    Ut_x1(imc,:) = phi(X1,:)*Utj;
%    Ut_x2(imc,:) = phi(X2,:)*Utj;
    
%    Ux_x1(imc,:) = grad_phi(X1,:)*Uj;
%    Ux_x2(imc,:) = grad_phi(X2,:)*Uj;
    
%    U_t1(imc,:) = phi*Utj(:,T1);
%    U_t2(imc,:) = phi*Utj(:,T2);
     
%    Ut_t1(imc,:) = phi*Utj(:,T1);
%    Ut_t2(imc,:) = phi*Utj(:,T2);
    
%    Ux_t1(imc,:) = grad_phi*Uj(:,T1);
%    Ux_t2(imc,:) = grad_phi*Uj(:,T2);

    % steady state vectors
      U_x1_ss = parfor_sub_vector( U_x1(imc,:),Tss);
%      U_x2_ss = parfor_sub_vector( U_x2(imc,:),Tss);
     Ut_x1_ss = parfor_sub_vector(Ut_x1(imc,:),Tss);
%     Ut_x2_ss = parfor_sub_vector(Ut_x2(imc,:),Tss);

    % compute FFT of u(X,t,w)
    U_x1_fft = ...
     sigproclib_fft_1sided(win'.*(U_x1_ss-mean(U_x1_ss)),freq_samp,nfft);
%    U_x2_fft = ...
%     sigproclib_fft_1sided(win'.*(U_x2_ss-mean(U_x2_ss)),freq_samp,nfft);

    % compute FFT of ut(X,t,w)
    Ut_x1_fft = ...
    sigproclib_fft_1sided(win'.*(Ut_x1_ss-mean(Ut_x1_ss)),freq_samp,nfft);
%    Ut_x2_fft = ...
%    sigproclib_fft_1sided(win'.*(Ut_x2_ss-mean(Ut_x2_ss)),freq_samp,nfft);
    
    % cumulative sum of FFT of u(X,t,w) square
    PSD_U_x1 = PSD_U_x1 + dt*abs(U_x1_fft).^2;
%    PSD_U_x2 = PSD_U_x2 + dt*abs(U_x2_fft).^2;
    
    % cumulative sum of FFT of ut(X,t,w) square
    PSD_Ut_x1 = PSD_Ut_x1 + dt*abs(Ut_x1_fft).^2;
%    PSD_Ut_x2 = PSD_Ut_x2 + dt*abs(Ut_x2_fft).^2;

    % L2 norm of the 2-norm of the unknown vector u(t)
    normL2_Q(1,imc) = norm_L2_time_vec(t0,t1,Uj);

    parfor_progress;
    
end
parfor_progress(0);


PSD_Ut_x1 = PSD_Ut_x1/Ns;
%PSD_Ut_x2 = PSD_Ut_x2/Ns;
%PSD_Ut_x3 = PSD_Ut_x3/Ns;



% compute MC convergence metric
conv_MC = MC_conv(normL2_Q,Ns);

clear Uj Utj;
clear U_x1_fft Ut_x1_fft;
clear U_x1_ss Ut_x1_ss;


% close MATLAB parallel pool
if parflag == 1
    disp(' ');
    matlabpool close
    disp(' ');
end


toc
% -----------------------------------------------------------



% save simulation data into a file
% -----------------------------------------------------------
tic

disp(' ')
disp(' --- Saving Workspace --- ');
disp(' ');
disp('    ... ');
disp(' ');


save([num2str(case_name),'.mat']);

toc
% -----------------------------------------------------------



% compute statistics
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- Computing Statistics --- ');
disp(' ');
disp('    ... ');
disp(' ');

aux_stat = {X1 X2 T1 T2 Tss ptile dt Ns Ndt nfft freq_samp win};
data_MC  = {U_x1 Ut_x1 ...
            PSD_U_x1 PSD_Ut_x1};

data_MC_stat = ...
randbar_fixed_mass_spring__statistics(aux_stat,data_MC);

clear aux_stat data_MC

toc
% -----------------------------------------------------------


% post processing
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- Post Processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


aux_post = {case_name t0 t1 Ndt Nx nfft sigma ptile Ns ...
            time xmesh freq phi WN Tss conv_MC};

randbar_fixed_mass_spring__post(aux_post,data_MC_stat);


toc
% -----------------------------------------------------------

