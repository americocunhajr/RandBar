
% ----------------------------------------------------------------- 
%  randbar_fixed_mass_spring__post.m
%
%  This script post process data that comes from the simulation 
%  of the stochastic dynamics of fixed-mass-spring bar.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Apr 6, 2013
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function randbar_fixed_mass_spring__post(aux_post,data_MC_stat)


% define auxiliar parameters
case_name = aux_post{1};
t0        = aux_post{2};
t1        = aux_post{3};
Ndt       = aux_post{4};
Nx        = aux_post{5};
nfft      = aux_post{6};
sigma     = aux_post{7};
ptile     = aux_post{8};
Ns        = aux_post{9};
time      = aux_post{10};
xmesh     = aux_post{11};
freq      = aux_post{12};
phi       = aux_post{13};
WN        = aux_post{14};
Tss       = aux_post{15};
conv_MC   = aux_post{16};


% define MC data statistics
% U_x1_mean      = data_MC_stat{1};
% U_x2_mean      = data_MC_stat{2};
%Ut_x1_mean      = data_MC_stat{3};
%Ut_x2_mean      = data_MC_stat{4};
%U_x1_upp        = data_MC_stat{5};
%U_x1_low        = data_MC_stat{6};
%U_x2_upp        = data_MC_stat{7};
%U_x2_low        = data_MC_stat{8};
%Ut_x1_upp       = data_MC_stat{9};
%Ut_x1_low       = data_MC_stat{10};
%Ut_x2_upp       = data_MC_stat{11};
%Ut_x2_low       = data_MC_stat{12};
%U_x1t1_pdf      = data_MC_stat{13};
%U_x1t1_supp     = data_MC_stat{14};
%U_x2t1_pdf      = data_MC_stat{15};
%U_x2t1_supp     = data_MC_stat{16};
%Ut_x1t1_pdf     = data_MC_stat{17};
%Ut_x1t1_supp    = data_MC_stat{18};
%Ut_x2t1_pdf     = data_MC_stat{19};
%Ut_x2t1_supp    = data_MC_stat{20};
%PSD_U_x1          = data_MC_stat{21};
%PSD_U_x2          = data_MC_stat{22};
%PSD_Ut_x1         = data_MC_stat{23};
%PSD_Ut_x2         = data_MC_stat{24};
%U_x1_nominal    = data_MC_stat{25};
%U_x2_nominal    = data_MC_stat{26};
%Ut_x1_nominal   = data_MC_stat{27};
%Ut_x2_nominal   = data_MC_stat{28};
%PSD_U_x1_nominal  = data_MC_stat{29};
%PSD_U_x2_nominal  = data_MC_stat{30};
%PSD_Ut_x1_nominal = data_MC_stat{31};
%PSD_Ut_x2_nominal = data_MC_stat{32};

% define MC data statistics
 U_x1_mean        = data_MC_stat{1};
Ut_x1_mean        = data_MC_stat{2};
U_x1_upp          = data_MC_stat{3};
U_x1_low          = data_MC_stat{4};
Ut_x1_upp         = data_MC_stat{5};
Ut_x1_low         = data_MC_stat{6};
U_x1t1_pdf        = data_MC_stat{7};
U_x1t1_supp       = data_MC_stat{8};
Ut_x1t1_pdf       = data_MC_stat{9};
Ut_x1t1_supp      = data_MC_stat{10};
PSD_U_x1          = data_MC_stat{11};
PSD_Ut_x1         = data_MC_stat{12};
U_x1_nominal      = data_MC_stat{13};
Ut_x1_nominal     = data_MC_stat{14};
PSD_U_x1_nominal  = data_MC_stat{15};
PSD_Ut_x1_nominal = data_MC_stat{16};
U_x1t1_bins       = data_MC_stat{17};
U_x1t1_freq       = data_MC_stat{18};
Ut_x1t1_bins      = data_MC_stat{19};
Ut_x1t1_freq      = data_MC_stat{20};



U_x1_nominal_ss = parfor_sub_vector(U_x1_nominal,Tss);
%U_x2_nominal_ss = parfor_sub_vector(U_x2_nominal,Tss);
%U_x3_nominal_ss = parfor_sub_vector(U_x3_nominal,Tss);

Ut_x1_nominal_ss = parfor_sub_vector(Ut_x1_nominal,Tss);
%Ut_x2_nominal_ss = parfor_sub_vector(Ut_x2_nominal,Tss);
%Ut_x3_nominal_ss = parfor_sub_vector(Ut_x3_nominal,Tss);


U_x1_mean_ss = parfor_sub_vector(U_x1_mean,Tss);
%U_x2_mean_ss = parfor_sub_vector(U_x2_mean,Tss);
%U_x3_mean_ss = parfor_sub_vector(U_x3_mean,Tss);

Ut_x1_mean_ss = parfor_sub_vector(Ut_x1_mean,Tss);
%Ut_x2_mean_ss = parfor_sub_vector(Ut_x2_mean,Tss);
%Ut_x3_mean_ss = parfor_sub_vector(Ut_x3_mean,Tss);




% plot U(L,t1,w) PDF
% -----------------------------------------------------------
gtitle = [' lumped mass displacement at t= ',num2str(t1,'%1.3f'),' sec'];
xlab   = ' normalized displacement';
ylab   = ' probability density function';
xmin   = -4.0;
xmax   =  4.0;
ymin   = 0.0;
ymax   = 2.5;
gname  = [num2str(case_name),'__pdf_u_x1t1'];
flag   = 'eps';
%fig1a  = graph_type1(U_x1t1_supp,U_x1t1_pdf,gtitle,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1a  = graph_bar_curve1(U_x1t1_bins,U_x1t1_freq,...
                          U_x1t1_supp,U_x1t1_pdf,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig1a);
% -----------------------------------------------------------


% plot U(X2,T1,w) PDF
% -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' norm. displac. at $x=2L/3$ and $t=T$';
% ylab   = ' probability density function';
% xmin   = -4.0;
% xmax   =  4.0;
% ymin   = 0.0;
% ymax   = 1.2;
% gname  = [num2str(case_name),'__pdf_u_x2t1'];
% flag   = 'eps';
%fig1b = graph_type1(U_x2t1_supp,U_x2t1_pdf,gtitle,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1b);
% -----------------------------------------------------------


% plot Ut(L,t1,w) PDF
% -----------------------------------------------------------
gtitle = [' lumped mass velocity at t= ',num2str(t1,'%1.3f'),' sec'];
xlab   = ' normalized velocity';
ylab   = ' probability density function';
xmin   = -4.0;
xmax   =  4.0;
ymin   = 0.0;
ymax   = 4.0;
gname  = [num2str(case_name),'__pdf_ut_x1t1'];
flag   = 'eps';
%fig1d  = graph_type1(Ut_x1t1_supp,Ut_x1t1_pdf,gtitle,...
%                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
fig1d  = graph_bar_curve1(Ut_x1t1_bins,Ut_x1t1_freq,...
                          Ut_x1t1_supp,Ut_x1t1_pdf,gtitle,...
                          xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig1d);
% -----------------------------------------------------------


% plot Ut(X2,T1,w) PDF
% -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' norm. veloc. at $x=2L/3$ and $t=T$';
% ylab   = ' probability density function';
% xmin   = -6.0;
% xmax   =  6.0;
% ymin   = 0.0;
% ymax   = 1.2;
% gname  = [num2str(case_name),'__pdf_ut_x2t1'];
% flag   = 'eps';
%fig1e = graph_type1(Ut_x2t1_supp,Ut_x2t1_pdf,gtitle,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig1e);
% -----------------------------------------------------------


% plot kNL(w) PDF
% -----------------------------------------------------------
% gtitle = ' ';
% xlab   = ' nonlinear stifeness $(N/m^3)$';
% ylab   = ' probability density function';
% xmin   = 'auto';
% xmax   = 'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% gname  = [num2str(case_name),'__pdf_kNL'];
% flag   = 'eps';
% fig2   = graph_type1(kNL_supp,kNL_pdf,gtitle,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
% close(fig2);
% -----------------------------------------------------------


% plot U(L,t,w) confidence interval
% -----------------------------------------------------------
gtitle = ' lumped mass displacement';
leg1   = ' deterministic';
leg2   = ' mean value';
leg3   = [' ',num2str(2*ptile-100),'% prob.'];
xlab   = ' time (ms)';
ylab   = ' displacement (\mum)';
xmin   = t0*1.0e3;
xmax   = t1*1.0e3;
ymin   = -600;
ymax   =  600;
gname  = [num2str(case_name),'__ci_u_x1'];
flag   = 'eps';
fig3a  = graph_ci2(time*1.0e3,U_x1_nominal*1.0e6,...
                   time*1.0e3,U_x1_mean*1.0e6,...
                   U_x1_upp*1.0e6,U_x1_low*1.0e6,gtitle,leg1,leg2,leg3,...
                   xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig3a);
% -----------------------------------------------------------


% plot Ut(L,t,w) confidence interval
% -----------------------------------------------------------
gtitle = ' lumped mass velocity';
leg1   = ' deterministic';
leg2   = ' mean value';
leg3   = [' ',num2str(2*ptile-100),'% prob.'];
xlab   = ' time (ms)';
ylab   = ' velocity (m/s)';
xmin   = t0*1.0e3;
xmax   = t1*1.0e3;
ymin   = -15;
ymax   =  15;
gname  = [num2str(case_name),'__ci_ut_x1'];
flag   = 'eps';
fig3b  = graph_ci2(time*1.0e3,Ut_x1_nominal,...
                   time*1.0e3,Ut_x1_mean,...
                   Ut_x1_upp,Ut_x1_low,gtitle,leg1,leg2,leg3,...
                   xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig3b);
% -----------------------------------------------------------



% plot U(X2,t,w) confidence interval
% -----------------------------------------------------------
% gtitle = ' ';
% leg1  = ' deterministic';
% leg2  = ' mean value';
% leg3  = [' ',num2str(2*ptile-100),'\%~prob.~~~'];
% xlab  = ' time (ms)';
% ylab  = ' displac. at $x=2L/3~~(\mu m)$';
% xmin  = t0*1.0e3;
% xmax  = t1*1.0e3;
% ymin   = 'auto';
% ymax   = 'auto';
% gname = [num2str(case_name),'__ci_u_x2'];
% flag  = 'eps';
%fig3c  = graph_ci2(time*1.0e3,U_x2_nominal*1.0e6,...
%                   time*1.0e3,U_x2_mean*1.0e6,...
%                   U_x2_upp*1.0e6,U_x2_low*1.0e6,gtitle,leg1,leg2,leg3,...
%                   xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3c);
% -----------------------------------------------------------


% plot Ut(X2,t,w) confidence interval
% -----------------------------------------------------------
% gtitle = ' ';
% leg1  = ' deterministic';
% leg2  = ' mean value';
% leg3  = [' ',num2str(2*ptile-100),'\% prob.'];
% xlab  = ' time~~$(ms)$';
% ylab  = ' veloc. at $x=2L/3~~(m/s)$';
% xmin  = t0*1.0e3;
% xmax  = t1*1.0e3;
% ymin   = 'auto';
% ymax   = 'auto';
% gname = [num2str(case_name),'__ci_ut_x2'];
% flag  = 'eps';
%fig3d = graph_ci2(time*1.0e3,Ut_x2_nominal,...
%                  time*1.0e3,Ut_x2_mean,...
%                   Ut_x2_upp,Ut_x2_low,gtitle,leg1,leg2,leg3,...
%                   xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig3d);
% -----------------------------------------------------------


% plot phase space at x=L
% -----------------------------------------------------------
gtitle = ' lumped mass orbit in phase space';
leg1   = ' deterministic';
leg2   = ' mean value';
xlab   = ' displacement (\mum)';
ylab   = ' velocity (m/s)';
xmin   = -600;
xmax   =  600;
ymin   = -15;
ymax   =  15;
gname  = [num2str(case_name),'__phase_space_x1'];
flag   = 'eps';
fig4a  = graph_type2(U_x1_nominal_ss*1.0e6,Ut_x1_nominal_ss,...
                     U_x1_mean_ss*1.0e6,Ut_x1_mean_ss,...
                     gtitle,leg1,leg2,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig4a);
% -----------------------------------------------------------


% plot phase space at x=X2
% -----------------------------------------------------------
% gtitle = ' ';
% leg1  = ' deterministic';
% leg2  = ' mean value';
% xlab  = ' displac. at $x=2L/3~~(\mu m)$';
% ylab  = ' veloc. at $x=2L/3~~(m/s)$';
% xmin  =  'auto';
% xmax  =  'auto';
% ymin   = 'auto';
% ymax   = 'auto';
% gname = [num2str(case_name),'__phase_space_x2'];
% flag  = 'eps';
%fig4b = graph_type2(U_x2_nominal_ss*1.0e6,Ut_x2_nominal_ss,...
%                    U_x2_mean_ss*1.0e6,Ut_x2_mean_ss,...
%                    gtitle,leg1,leg2,...
%                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig4b);
% -----------------------------------------------------------


% external force at x=L
% -----------------------------------------------------------
gtitle = ' external force at the lumped mass';
xlab  = ' time (ms)';
ylab  = ' external force (kN)';
xmin  = t0*1.0e3;
xmax  = t1*1.0e3;
ymin  = 'auto';
ymax  = 'auto';
gname = [num2str(case_name),'__F_x1'];
flag  = 'eps';
fig5a = graph_typeN(time*1.0e3,sigma*phi(Nx,1).*WN(1,:)*1.0e-3,...
                  gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

close(fig5a);
% -----------------------------------------------------------


% external force at t=T
% -----------------------------------------------------------
% gtitle = ' ';
% xlab  = ' position$(m)$';
% ylab  = ' external force at $t=T (kN)$';
% xmin  = 0.0;
% xmax  = 1.0;
% ymin  = 'auto';
% ymax  = 'auto';
% gname = [num2str(case_name),'__F_T1'];
% flag  = 'eps';
%fig5b = graph_type1(xmesh,sigma*phi(:,1).*WN(Ns,Ndt)*1.0e-3,...
%                   gtitle,xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);

%close(fig5b);
% -----------------------------------------------------------



% displacement spectral density at x=L
% -----------------------------------------------------------
gtitle = ' lumped mass displacement';
leg1  = ' deterministic';
leg2  = ' mean value';
xlab  = ' frequency (kHz)';
ylab  = ' power spectral density (dB/kHz)';
xmin  = 0.0;
xmax   = 'auto';
ymin   = 'auto';
ymax   = 'auto';
gname = [num2str(case_name),'__PSD_u_x1'];
flag  = 'eps';
fig6a = graph_type2(freq*1.0e-3,10*log10(PSD_U_x1_nominal*1.0e3),...
                    freq*1.0e-3,10*log10(PSD_U_x1*1.0e3),...
                    gtitle,leg1,leg2,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig6a);
% -----------------------------------------------------------


% velocity spectral density at x=L
% -----------------------------------------------------------
gtitle = ' lumped mass velocity';
leg1   = ' deterministic';
leg2   = ' mean value';
xlab   = ' frequency (kHz)';
ylab   = ' power spectral density (dB/kHz)';
xmin   = 0.0;
ymin   = 'auto';
ymax   = 'auto';
ymax   = 0.0;
gname  = [num2str(case_name),'__PSD_ut_x1'];
flag   = 'eps';
fig6b  = graph_type2(freq*1.0e-3,10*log10(PSD_Ut_x1_nominal*1.0e3),...
                     freq*1.0e-3,10*log10(PSD_Ut_x1*1.0e3),...
                     gtitle,leg1,leg2,...
                     xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig6b);
% -----------------------------------------------------------


% displacement spectral density at x=X2
% -----------------------------------------------------------
% gtitle = ' ';
% leg1  = ' deterministic';
% leg2  = ' mean value';
% xlab  = ' frequency (Hz)';
% ylab  = ' displac. PSD at $x=2L/3$';
% xmin  = 'auto';
% xmax  = 'auto';
% ymin  = 'auto';
% ymax  = 'auto';
% gname = [num2str(case_name),'__PSD_u_x2'];
% flag  = 'eps';
%fig6c = graph_loglog2(freq,PSD_U_x2,freq,PSD_U_x2_nominal,gtitle,leg1,leg2,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6c);
% -----------------------------------------------------------


% velocity spectral density at x=X2
% -----------------------------------------------------------
% gtitle = ' ';
% leg1  = ' deterministic';
% leg2  = ' mean value';
% xlab  = ' frequency (Hz)';
% ylab  = ' veloc. PSD at $x=L$';
% xmin  = 'auto';
% xmax  = 'auto';
% ymin  = 'auto';
% ymax  = 'auto';
% gname = [num2str(case_name),'__PSD_ut_x2'];
% flag  = 'eps';
%fig6d = graph_loglog2(freq,PSD_Ut_x2,freq,PSD_Ut_x2_nominal,gtitle,leg1,leg2,...
%                      xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
%close(fig6d);
% -----------------------------------------------------------


% MC convergence metric
% -----------------------------------------------------------
gtitle = ' MC convergence metric';
xlab   = ' number of MC realizations';
ylab   = ' convergence metric';
xmin   = 0;
xmax   = Ns;
ymin   = 'auto';
ymax   = 'auto';
gname  = [num2str(case_name),'__MC_conv'];
flag   = 'eps';
fig7   = graph_type1((1:Ns),conv_MC,gtitle,...
                    xlab,ylab,xmin,xmax,ymin,ymax,gname,flag);
close(fig7);
% -----------------------------------------------------------


% probability of U(L,t1,.) <= 0
U_x1t1_supp_int = abs(U_x1t1_supp.*(U_x1t1_supp <= 0.0));
[value index] = min(U_x1t1_supp_int);
prob_U_x1t1 = trapz(U_x1t1_supp(1:index),U_x1t1_pdf(1:index))



% probability of Ut(L,t1,.) <= 0
Ut_x1t1_supp_int = abs(Ut_x1t1_supp.*(Ut_x1t1_supp <= 0.0));
[value index] = min(Ut_x1t1_supp_int);
prob_Ut_x1t1 = trapz(Ut_x1t1_supp(1:index),Ut_x1t1_pdf(1:index))


% skewness of U(L,t1,.)
%skewness()


return
% -----------------------------------------------------------------

