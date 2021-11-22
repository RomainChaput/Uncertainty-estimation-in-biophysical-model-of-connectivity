%% Code to Analyze the outputs of the Sparse Pseudo-Spectral Projection - 2020/04/27
% Adapted from M. Iskandarani and Pierre Sochala
%
% Author: R. Chaput - Modified for publication on 2021/11/22


clear
close

%% First initialize the PC-information and form the projection matrix for level 3

addpath(genpath('./PSP'));
PC_data = PSP_structure(5,3); % 5 inputs and 3 levels

% load the input values and order of runs
load('Level3_input_coord.mat');

%Input_values = [Settlement_speed_input; Kappa_input; Competency_input; Flexion_input; Beta_input]';
pmin=[min(Settlement);min(Kappa);min(Competency);min(Flexion)/(3600*24);min(Beta)];   % minima of uncertain parameters
pmax=[max(Settlement);max(Kappa);max(Competency);max(Flexion)/(3600*24);max(Beta)];   % maxima of uncertain parameters

% Read-in the realizations
addpath(genpath('./Dispersal_kernels_samples'));
load('dispersal_kernel_distance_lower_keys.dat');
% Results for the 4 reefs of interest
fr_lower_keys = (dispersal_kernel_distance_lower_keys(:,2:352)); % First column gives the spatial points of the dispersal kernels

xp = (0:5:350);% location of spatial points in the dispersal kernels
nup = 95; % number of polynomials
psinor = PC_data.sqnorm; % weighed 2-norm of multi-dimensional basis
Anisp = PC_data.Matrix; % projection matrix of realization
pcnpt = PC_data.Multi_ind; % multi-index of polynomial basis                 

%% Compute PC coefficients using Least Square

% Least square method is prefered here to the projection method to remove
% output noise
N_SP = 351; % Number of sample CMS runs 
SP = X5D_SG_level3(1:N_SP,:)*2-1; % Sparse grid points
M = evaluation_PC(SP,pcnpt);
A_ls = (M'*M)\(M');
% PC coefficients
fh_lower = (fr_lower_keys*A_ls'); % produce the PC coefficients \hat{f}_n
                                  % fh(:,1) is the 0-th mode or the mean
                                  % fh(:,m) is the (m-1)-th multi-dimensional mode, => mean of the output
                                  % m=1,...nup+1 
                                  
% % Compute PC-coefficients using projection method
% fh_lower = (fr_lower_keys*Anisp');                 

% Compute the standard deviation of the modal output at all the points of interest, for all the variables
fstdev_lower = sqrt( (fh_lower(:,2:nup+1).^2) * psinor(2:nup+1) );


%% Use the PC surrogates to plot mean dispersal kernel and plus or minus a standard deviation

% Plot mean and plus or minus a standard deviation
addpath(genpath('./Output_figures'));
figure()
plot(xp, fh_lower(:,1),'k',... % First column of the surrogate is the mean PC dispersal kernel
     xp, fh_lower(:,1)+fstdev_lower,'r--',...
     xp, fh_lower(:,1)-fstdev_lower,'r--');
 hold on
 plot(xp, dispersal_kernel_distance_lower_keys(:,2),'b--'); % Dispersal kernel as estimated by the CMS with mean input values
 hl=legend({'$\hat{M}$','$\hat{M}+\sigma(x)$','$\hat{M}-\sigma(x)$','CR'},'Interpreter','latex','location','northeast');
 set(hl,'FontSize',16);
 ylim([0 8*10^(-4)])
 xlabel('Distance (Km)','FontSize',12);
 ylabel('Proportion of settlers per reef polygon','fontsize',12);
 title('PC Dispersal Kernels - Lower Keys')
savefig('./Output_figures/Lower Keys dispersal kernel.fig')


%% Validation points for the PC surrogate

% Load realizations from CMS simulations for validation points
addpath(genpath('./Validation_dispersal_kernels'));
N_VP = 100; % Number of validation points

% Results of the CMS runs with validation input parameters for Lower Keys
load('validation_dispersal_kernel_distance_lower_keys.dat') 
DK_VP_CMS = (validation_dispersal_kernel_distance_lower_keys(:,2:(N_VP+1))); % CMS dispersal kernels at the validation points

% Computation of PC estimates for the validation points
fh = fh_lower; % PC coefficients for the Lower Keys (built using the 351 sample points)
load('Validation_points_input_coord.mat') % Values of input parameters for validation runs
VP = Validation_points_coord(1:100,:)*2-1; % Normalized validation point coordinates
clear Validation_points_coord
ndim = 5; % Number of dimensions
nord = 7; % Order of the polynomials
DK_VP_PC = zeros(length(fh(:,1)),N_VP);

for iv = 1:N_VP
    Pv1d1 = zeros(nord+1,nord+1);
    for j = 1:5
        Pv1d1(:,j) = legendrepols(nord,VP(iv,j)); % Computation of 1D polynomials for 1 validation point for each of the 5 coordinates
    end
    Pvnd = ones([nup+1 1]);
    for p = 1:nup
        for d = 1:ndim
            n = pcnpt(p,d);
            Pvnd(p) = Pvnd(p) * Pv1d1(n+1,d); % Multiplicative loop over the legendre polynomes for the 5 dimensions
        end
    end
    DK_VP_PC(:,iv) = fh * Pvnd; %PC estimate for the dispersal kernels at the validation points
end


%% Computation of the fit of the model: RMSE

RMSE = sqrt(mean((DK_VP_PC' - DK_VP_CMS').^2));
figure()
plot(xp,RMSE);
ylabel('RMSE')
xlabel('Dispersal Distance (km)')
title('RMSE - Lower Keys')%
savefig('./Output_figures/LK_RMSE_5D_100VP.fig')
%
figure()
RMSE_normalized = RMSE./fh(:,1)';
%RMSE_normalized = RMSE./(mean(DK_VP_CMS'));
plot(xp,RMSE_normalized);
xlim([0 350])
ylabel('Scatter Index (Normalized RMSE)')
xlabel('Dispersal Distance (km)')
title('RMSE Normalized - Lower Keys')%
savefig('./Output_figures/LK_RMSE_normalized_5D_100VP.fig')


%% Analysis of the Residuals as function of the input parameters

Residuals = (DK_VP_PC - DK_VP_CMS)';
xlabs = {'Sw.S.', 'Kappa', 'Comp.', 'Flex.', 'Beta'}; % Names of the uncertain input parameters
figure()
for i = 1:5
    subplot(2,5,i)
    scatter(VP(:,i),mean(Residuals,2),'.')
    hold on 
    p=polyfit(VP(:,i),mean(Residuals,2),1);
    f = polyval(p,VP(:,i)); 
    plot(VP(:,i),f,'-')
    [R,P] = corrcoef(VP(:,i),mean(Residuals,2));
    Rsq = R(1,2).^2;
    plot(0, 0, '.k','MarkerSize',1)
    hold on
    legend(sprintf('Residuals'),sprintf('Linear Regression'),sprintf('R^2 = %.2f',Rsq))
    xlabel(xlabs{i})
    hline(0,'k');
end
subplot(2,5,1)
ylabel('Mean Residuals')
for i = 6:10
    subplot(2,5,i)
    scatter(VP(:,i-5),std(Residuals,0,2)','.')
    hold on
    p=polyfit(VP(:,i-5),std(Residuals,0,2)',1);
    f = polyval(p,VP(:,i-5)); 
    plot(VP(:,i-5),f,'-')
    [R,P] = corrcoef(VP(:,i-5),std(Residuals,0,2)');
    Rsq = R(1,2).^2;
    plot(0, 0, '.k','MarkerSize',1)
    hold on
    legend(sprintf('Residuals'),sprintf('Linear Regression'),sprintf('R^2 = %.2f',Rsq))
    xlabel(xlabs{i-5})
    hline(0,'k');
end
subplot(2,5,6)
ylabel('Standard Deviation Residuals')
subplot(2,5,8)
xlabel({'Comp.';'Standardized Input Parameters'})
savefig('./Output_figures/LK_residual_dist.fig')

% Quantile plots: Residuals VS Normal distribution and surrogate VS CMS
figure()
pd = makedist('Normal');
qqplot(mean(Residuals,2),pd)
title('QQplot - Lower Keys')%
savefig('./Output_figures/LK_qqplot_Residuals_Normal_dist.fig')
%
figure()
qqplot(mean(DK_VP_PC),mean(DK_VP_CMS))
title('QQplot PC/CMS - Lower Keys')%
savefig('./Output_figures/LK_qqplot_PC_CMS.fig')

%% PC coefficients plots

% Drawing a 3D bar plot to show coefficients
figure()
ip = 1;% distance for analysis: 1,3,11,21 (0, 10, 50, 100km)
fh2d = nan([nord+1 nord+1]);
for i=1:nup+1
 id = pcnpt(i,1);   % polynomial degree along dimension 1
 jd = pcnpt(i,2);   % polynomial degree along dimension 2
 fh2d(id+1,jd+1) = fh_lower(ip,i);
end
xbar=0:1:nord;
bar3(abs(fh2d)');
xlabel('$m+1$','interpret','latex');
ylabel('$n+1$','interpret','latex');
zlabel('$|\widehat{M}_{m,n}|$','interpret','latex');
set(gca,'FontName','Times','FontSize',16);
savefig('./Output_figures/Lower Keys spectral coefficients at 0km.fig')


%% Compute the Sobol indices to analyze the variability induced by the input parameters

% Squared total variance
fh = fh_lower;
variance2 = (fh(:,2:96)).^2;
D = sum(variance2,2);

% Computation of First order Sobol indices: Si
Si = zeros(5,71);
for i = 1:5  
    pcnpt_index = pcnpt;
    pcnpt_index(:,i)=[];
    idx_i = find(all(pcnpt_index'==0)); 
    for dist = 1:71
        Si(i,dist) = (sum(fh(dist,idx_i(2:5)).^2))/D(dist);
    end
end

% Total Sobol indices per input parameter
Ti = zeros(5,71);
for i = 1:5
    pcnpt_index = pcnpt(:,i);
    idx_total_order = find(pcnpt_index~=0);
    for dist = 1:71
        Ti(i,dist) = (sum(fh(dist,idx_total_order(1:length(idx_total_order))).^2))/D(dist);
    end
end
% Computation of the fraction of variability introduced by each parameter
sumTi = sum(Ti); 
fracTi = Ti./sumTi;

% Short distance of dispersal
shortTi2=zeros(1,5);
for i = 1:5
    shortTi2(i) = mean(rmmissing(fracTi(i,1:3)));
end
% Medium distance of dispersal
medTi2=zeros(1,5);
for i = 1:5
    medTi2(i) = mean(rmmissing(fracTi(i,4:21)));
end
% Long distance of dispersal
longTi2=zeros(1,5);
for i = 1:5
    longTi2(i) = mean(rmmissing(fracTi(i,21:71)));
end

% Si as fraction of Ti
fracSi = Si./sumTi;
shortSi2=zeros(1,5);
for i = 1:5
    shortSi2(i) = mean(rmmissing(fracSi(i,1:3)));
end
% Medium distance of dispersal
medSi2=zeros(1,5);
for i = 1:5
    medSi2(i) = mean(rmmissing(fracSi(i,4:21)));
end
% Long distance of dispersal
longSi2=zeros(1,5);
for i = 1:5
    longSi2(i) = mean(rmmissing(fracSi(i,21:71)));
end

% Lower Keys
subplot(1,2,1)
x = categorical({'0-10km','10-100km','>100km'});
vals = [shortSi2;medSi2;longSi2];
bar(x,vals)
ylabel('Fraction of variance')
title('A. First order Lower Keys')
ylim([0 1])
legend('Sw.S.','Kappa','Comp.','Flex.','Beta')
%
subplot(1,2,2)
x = categorical({'0-10km','10-100km','>100km'});
vals = [shortTi2;medTi2;longTi2];
bar(x,vals)
title('B. Total order Lower Keys')
ylim([0 1])
legend('Sw.S.','Kappa','Comp.','Flex.','Beta')
yyaxis right
set(gca,'YTickLabel',[])
set(gca,'YTick')
set(gca,'ycolor','k') 
savefig('./Output_figures/Sobol_indices_LK.fig')