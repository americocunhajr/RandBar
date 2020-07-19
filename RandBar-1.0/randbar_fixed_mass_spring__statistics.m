
% ----------------------------------------------------------------- 
%  randbar_fixed_mass_spring__statistics.m
%
%  This function computes statistics of the data that comes from 
%  the simulation of fixed-mass-spring bar stochastic dynamics.
% ----------------------------------------------------------------- 
%  programmer: Americo Barbosa da Cunha Junior
%              americo.cunhajr@gmail.com
%
%  last update: Apr 6, 2013
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function data_MC_stat = ...
         randbar_fixed_mass_spring__statistics(aux_stat,data_MC)


% define auxiliar parameters
X1           = aux_stat{1};
X2           = aux_stat{2};
T1           = aux_stat{3};
T2           = aux_stat{4};
Tss          = aux_stat{5};
ptile        = aux_stat{6};
dt           = aux_stat{7};
Ns           = aux_stat{8};
Ndt          = aux_stat{9};
nfft         = aux_stat{10};
freq_samp    = aux_stat{11};
win          = aux_stat{12};

% define MC data matrices
% U_x1   = data_MC{1};
% U_x2   = data_MC{2};
%Ut_x1   = data_MC{3};
%Ut_x2   = data_MC{4};
%PSD_U_x1  = data_MC{5};
%PSD_U_x2  = data_MC{6};
%PSD_Ut_x1 = data_MC{7};
%PSD_Ut_x2 = data_MC{8};

 U_x1     = data_MC{1};
Ut_x1     = data_MC{2};
PSD_U_x1  = data_MC{3};
PSD_Ut_x1 = data_MC{4};


% nominal values
U_x1_nominal = U_x1(Ns,:);
%U_x2_nominal = U_x2(Ns,:);

Ut_x1_nominal = Ut_x1(Ns,:);
%Ut_x2_nominal = Ut_x2(Ns,:);

U_x1_nominal_ss = parfor_sub_vector(U_x1_nominal(1,:),Tss);
%U_x2_nominal_ss = parfor_sub_vector(U_x2_nominal(1,:),Tss);

Ut_x1_nominal_ss = parfor_sub_vector(Ut_x1_nominal(1,:),Tss);
%Ut_x2_nominal_ss = parfor_sub_vector(Ut_x2_nominal(1,:),Tss);

U_x1_nominal_fft = ...
     sigproclib_fft_1sided(win'.*(U_x1_nominal_ss-mean(U_x1_nominal_ss)),freq_samp,nfft);
%U_x2_nominal_fft = ...
%     sigproclib_fft_1sided(win'.*(U_x2_nominal_ss-mean(U_x2_nominal_ss)),freq_samp,nfft);

Ut_x1_nominal_fft = ...
     sigproclib_fft_1sided(win'.*(Ut_x1_nominal_ss-mean(Ut_x1_nominal_ss)),freq_samp,nfft);
%Ut_x2_nominal_fft = ...
%     sigproclib_fft_1sided(win'.*(Ut_x2_nominal_ss-mean(Ut_x2_nominal_ss)),freq_samp,nfft);

PSD_U_x1_nominal = dt*abs(U_x1_nominal_fft).^2;
%PSD_U_x2_nominal = dt*abs(U_x2_nominal_fft).^2;

PSD_Ut_x1_nominal = dt*abs(Ut_x1_nominal_fft).^2;
%PSD_Ut_x2_nominal = dt*abs(Ut_x2_nominal_fft).^2;



% estimate mean value
U_x1_mean = mean(U_x1);
%U_x2_mean = mean(U_x2);

Ut_x1_mean = mean(Ut_x1);
%Ut_x2_mean = mean(Ut_x2);


% estimate standard deviation
% U_x1_std = std(U_x1);
% U_x2_std = std(U_x2);

% Ut_x1_std = std(Ut_x1);
% Ut_x2_std = std(Ut_x2);



% estimate percentiles
U_x1_upp = prctile(U_x1,ptile);
%U_x2_upp = prctile(U_x2,ptile);
% U_x3_upp = prctile(U_x3,ptile);

U_x1_low = prctile(U_x1,100-ptile);
%U_x2_low = prctile(U_x2,100-ptile);
% U_x3_low = prctile(U_x3,100-ptile);


Ut_x1_upp = prctile(Ut_x1,ptile);
%Ut_x2_upp = prctile(Ut_x2,ptile);
%Ut_x3_upp = prctile(Ut_x3,ptile);

Ut_x1_low = prctile(Ut_x1,100-ptile);
%Ut_x2_low = prctile(Ut_x2,100-ptile);
% Ut_x3_low = prctile(Ut_x3,100-ptile);


% estimate probability density function
[U_x1t1_pdf, U_x1t1_supp] = ksdensity(sigproclib_norm_randvar(U_x1(:,T1)));
[U_x1t1_bins,U_x1t1_freq] = sigproclib_pdf1(sigproclib_norm_randvar(U_x1(:,T1)),sqrt(Ns),-4.0,4.0);
%[U_x2t1_pdf,U_x2t1_supp] = ksdensity(sigproclib_norm_randvar(U_x2(:,T1)));

%[U_x1t2_pdf,U_x1t2_supp] = ksdensity(sigproclib_norm_randvar(U_x1(:,T2)));
%[U_x2t2_pdf,U_x2t2_supp] = ksdensity(sigproclib_norm_randvar(U_x2(:,T2)));

[Ut_x1t1_pdf, Ut_x1t1_supp] = ksdensity(sigproclib_norm_randvar(Ut_x1(:,T1)));
[Ut_x1t1_bins,Ut_x1t1_freq] = sigproclib_pdf1(sigproclib_norm_randvar(Ut_x1(:,T1)),sqrt(Ns),-4.0,4.0);
%[Ut_x2t1_pdf,Ut_x2t1_supp] = ksdensity(sigproclib_norm_randvar(Ut_x2(:,T1)));

%[Ut_x1t2_pdf,Ut_x1t2_supp] = ksdensity(sigproclib_norm_randvar(Ut_x1(:,T2)));
%[Ut_x2t2_pdf,Ut_x2t2_supp] = ksdensity(sigproclib_norm_randvar(Ut_x2(:,T2)));

%[Ux_x1t1_pdf,Ux_x1t1_supp] = ksdensity(sigproclib_norm_randvar(Ux_x1(:,T1)));
%[Ux_x2t1_pdf,Ux_x2t1_supp] = ksdensity(sigproclib_norm_randvar(Ux_x2(:,T1)));

%[Ux_x1t2_pdf,Ux_x1t2_supp] = ksdensity(sigproclib_norm_randvar(Ux_x1(:,T2)));
%[Ux_x2t2_pdf,Ux_x2t2_supp] = ksdensity(sigproclib_norm_randvar(Ux_x2(:,T2)));

%[E_pdf,E_supp] = ksdensity(E(:,1));
%[kNL_pdf,kNL_supp] = ksdensity(kNL(:,1));


% estimate spectral density function
PSD_U_x1(2:end-1) = 2.0*PSD_U_x1(2:end-1);
%PSD_U_x2(2:end-1) = 2.0*PSD_U_x2(2:end-1);
PSD_Ut_x1(2:end-1) = 2.0*PSD_Ut_x1(2:end-1);
%PSD_Ut_x2(2:end-1) = 2.0*PSD_Ut_x2;

PSD_U_x1_nominal(2:end-1) = 2.0*PSD_U_x1_nominal(2:end-1);
%PSD_U_x2_nominal(2:end-1) = 2.0*PSD_U_x2_nominal(2:end-1);
PSD_Ut_x1_nominal(2:end-1) = 2.0*PSD_Ut_x1_nominal(2:end-1);
%PSD_Ut_x2_nominal(2:end-1) = 2.0*PSD_Ut_x2_nominal(2:end-1);


% data statistics
% data_MC_stat = ...
%  {U_x1_mean U_x2_mean Ut_x1_mean Ut_x2_mean ...
%   U_x1_upp U_x1_low U_x2_upp U_x2_low ...
%   Ut_x1_upp Ut_x1_low Ut_x2_upp Ut_x2_low ...
%   U_x1t1_pdf U_x1t1_supp U_x2t1_pdf U_x2t1_supp ...
%   Ut_x1t1_pdf Ut_x1t1_supp Ut_x2t1_pdf Ut_x2t1_supp ...
%   PSD_U_x1 PSD_U_x2 PSD_Ut_x1 PSD_Ut_x2,...
%   U_x1_nominal,U_x2_nominal,Ut_x1_nominal,Ut_x2_nominal,...
%   PSD_U_x1_nominal,PSD_U_x2_nominal,PSD_Ut_x1_nominal,PSD_Ut_x2_nominal};

data_MC_stat = ...
 {U_x1_mean Ut_x1_mean ...
  U_x1_upp U_x1_low ...
  Ut_x1_upp Ut_x1_low ...
  U_x1t1_pdf U_x1t1_supp ...
  Ut_x1t1_pdf Ut_x1t1_supp ...
  PSD_U_x1 PSD_Ut_x1 ,...
  U_x1_nominal,Ut_x1_nominal,...
  PSD_U_x1_nominal,PSD_Ut_x1_nominal,...
  U_x1t1_bins,U_x1t1_freq,...
  Ut_x1t1_bins,Ut_x1t1_freq};


return
% -----------------------------------------------------------------
