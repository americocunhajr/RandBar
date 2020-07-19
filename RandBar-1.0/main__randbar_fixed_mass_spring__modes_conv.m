
% ----------------------------------------------------------------- 
%  main__randbar_fixed_mass_spring__modes_conv.m
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
%  This script is the main file for a program that study
%  the convergence of solution approximation as function
%  of the number of modes used in the approximation.
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
disp(' Fixed-Mass-Spring Bar Modes Convergence            ')
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
case_name = 'modes_conv';

disp(' '); 
disp([' Case Name: ',num2str(case_name)]);
disp(' ');
% -----------------------------------------------------------



% pre processing
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Pre Processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% simulation parameters
[t0,t1,rho,E,d,A,L,c,k,kNL,m,...
 alpha1,alpha2,sigma,nlflag,...
 Nx,dx,xmesh,dt,time,Ndt,...
 Nmodes,wn,phi,grad_phi,lambda] = ...
          randbar_fixed_mass_spring__phys_param(case_name);

% external force (N)
Fwn = ones(1,Ndt);
        
% physical parameters vector
phys_param = [rho E A L k kNL alpha1 alpha2 sigma c m];


toc
% -----------------------------------------------------------



% Error Estimation
% -----------------------------------------------------------
tic
disp(' '); 
disp(' --- Modes Convergence Study --- ');
disp(' ');
disp('    ... ');
disp(' ');


% compute bins vector
Nmodes_bins = 1:1:Nmodes;

% dimension of bins vector
Nbins = length(Nmodes_bins);

% bins counter
ibin = 0;

% preallocate memory for L2 residual vector
ResL2 = zeros(Nbins,1);

% preallocate memory for H1 residual vector
ResH1 = zeros(Nbins,1);

% preallocate memory for U0 matrix
U0 = zeros(Nx,Ndt);

% preallocate memory for V0 matrix
Ux0 = zeros(Nx,Ndt);



% study of modes convergence
% -----------------------------------------------------------
for imodes=Nmodes_bins

    disp([' number of modes: ',num2str(imodes)]);
    
    
    % Galerkin matrices and vectors
    [M,C,K,F,u0,v0] = ...
       randbar_fixed_mass_spring__galerkin(imodes,time,...
                                             lambda,Fwn,phys_param);
    
    
    % numerical integration of Galerkin ODE system
    [Uj,~,~,~] = ...
        newmark_ode_solver(M,C,K,F,u0,v0,...
                dt,Ndt,Nx,imodes,phi(:,1:imodes),nlflag,phys_param);
    
    
    % compute U and Ux from coeficients and shape modes
    U  = phi(:,1:imodes)*Uj;
    %Ut = phi(:,1:imodes)*Utj;
    Ux = grad_phi(:,1:imodes)*Uj;
    
    
    % error estimation
    ibin = ibin + 1;
    
    ResL2(ibin,1) = norm_L2(0,L,U(:,Ndt)-U0(:,Ndt));
    
    ResH1(ibin,1) = norm_H1(0,L,U(:,Ndt)-U0(:,Ndt),Ux(:,Ndt)-Ux0(:,Ndt));
    
    U0  = U;
    Ux0 = Ux;

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



% post processing
% -----------------------------------------------------------
tic
disp(' ');
disp(' --- Post Processing --- ');
disp(' ');
disp('    ... ');
disp(' ');


% plot ResL2 vs Nmodes
% -----------------------------------------------------------
xlab  = ' number of modes';
ylab  = ' $L_2$ norm of residual';
xmin  = 'auto';
xmax  = 'auto';
ymin  = 'auto';
ymax  = 'auto';
gname  = [num2str(case_name),'__resL2_vs_Nmodes'];
flag  = 'eps';
fig1  = graph_type1x(Nmodes_bins,ResL2,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig1  = graph_semilogy1x(Nmodes_bins,ResL2,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig1);
% -----------------------------------------------------------


% plot ResH1 vs Nmodes
% -----------------------------------------------------------
xlab  = ' number of modes';
ylab  = ' $H^1$ norm of residual';
xmin  = 'auto';
xmax  = 'auto';
ymin  = 'auto';
ymax  = 'auto';
gname  = [num2str(case_name),'__resH1_vs_Nmodes'];
flag  = 'eps';
fig2  = graph_type1x(Nmodes_bins,ResH1,...
                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%fig2  = graph_semilogy1x(Nmodes_bins,ResH1,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig2);
% -----------------------------------------------------------


% modes vs residuals metrics table
% -----------------------------------------------------------
disp(' ')
disp(' Nmodes vs L2/H1 residuals')
disp(' ')

fid1 = fopen('residual.dat', 'w');
for ibin=1:1:Nbins

    disp([' ',num2str(Nmodes_bins(ibin), '% 03d'),...
          ' ',num2str(ResL2(ibin), '% 3.5e'),...
          ' ',num2str(ResH1(ibin), '% 3.5e')]);
    fprintf(fid1,'\n% 03d  % 3.5e  % 3.5e',...
            Nmodes_bins(ibin),ResL2(ibin),ResH1(ibin));
end
fclose(fid1);
% -----------------------------------------------------------


toc
% -----------------------------------------------------------
