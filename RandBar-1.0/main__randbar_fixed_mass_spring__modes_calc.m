
% ----------------------------------------------------------------- 
%  main__randbar_fixed_mass_spring__modes_calc.m
%
%                                 |--> u(x,t)
%                                  ___
% /|______________________________|   |      k     |\
% /|                              |   |---/\/\/\---|\
% /| rho, E, A                    | m |            |\
% /|______________________________|   |            |\
% /|                              |___|            |\
%
%  |----> x
%
%  rho A u_tt =  [E A u_x]_x + f        for 0 < x < L
%           u =  0                      for x = 0
%     E A u_x = -k.u - m.u_tt           for x = L
%           u =  u0                     for t = 0
%         u_t =  v0                     for t = 0
%
%  This script is the main file for a program that computes
%  the modes of fixed-mass-spring bar.
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
disp(' Fixed-Mass-Spring Bar Modes                        ')
disp('                                                    ')
disp('                                                    ')
disp('                                 |--> u(x,t)        ')
disp('                                  ___               ')
disp(' /|______________________________|   |      k     |\')
disp(' /|                              |   |---/\/\/\---|\')
disp(' /| rho, E, A                    | m |            |\')
disp(' /|______________________________|   |            |\')
disp(' /|                              |___|            |\')
disp('                                                    ')
disp('  |----> x                                          ')
disp('                                                    ')
disp('  rho A u_tt =  [E A u_x]_x + f        for 0 < x < L')
disp('           u =  0                      for x = 0    ')
disp('     E A u_x = -k.u - m.u_tt           for x = L    ')
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
case_name = 'case1a';
%case_name = 'case1b';
%case_name = 'case1c';
%case_name = 'case1d';


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
 Nmodes,freq_nat,phi,grad_phi,lambda] = ...
          randbar_fixed_mass_spring__phys_param(case_name);
                                    
toc
% -----------------------------------------------------------


% save some variables to a file
% -----------------------------------------------------------
tic

disp(' ')
disp(' Saving Workspace: start ')


save('bar_modes.mat');

disp(' ')
disp(' Saving Workspace: end ')
disp(' ')

toc
% -----------------------------------------------------------



% post processing
% -----------------------------------------------------------
tic

disp(' ')
disp(' Post Processing: start ')


for imode=1:1:Nmodes


% plot phi(x) vs x
% -----------------------------------------------------------
if imode == 1
    gtitle = [num2str(imode),'-st mode'];
elseif imode == 2
    gtitle = [num2str(imode),'-nd mode'];
else
    gtitle = [num2str(imode),'-th mode'];
end
xlab   = ' x (m)';
ylab   = [' \phi_{',num2str(imode),'}'];
xmin   = 0.0;
xmax   = L;
ymin   = 'auto';
ymax   = 'auto';
gname  = ['bar_phi_vs_x_',num2str(imode)];
flag   = 'eps';
fig1   = graph_type1(xmesh,phi(:,imode),gtitle,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig1);
% -----------------------------------------------------------

end


disp(' ')
disp(' Post Processing: end ')
disp(' ')

toc
% -----------------------------------------------------------
