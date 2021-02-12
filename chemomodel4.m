function chemomodel4

close all; clear all;

%
% MLMC treatment
%

addpath('..');

M     = 2;      % refinement cost factor
Lmin  = 3;      % minimum refinement level
Lmax  = 20;     % maximum refinement level

N      = 8000;    % samples for convergence tests
L      = 8;        % levels for convergence tests 
N0     = 200;      % initial samples on each level
%Eps    = [ 0.005 0.01 0.02 0.05];
%Eps    = [ 0.02 0.05 0.1 0.2 0.5];
%Eps    = [ 0.05 0.1 0.2 0.5 1];
Eps = [0.1 0.2 0.5 1];

filename = 'chemomodel4';
fp       = fopen([filename '.txt'],'w');
mlmc_test(@evsi_l, M, N,L, N0,Eps,Lmin,Lmax, fp);
fclose(fp);

%
% plot results
%

nvert = 2;
mlmc_plot(filename, nvert);

if(nvert==1)
  figure(1)
  print('-deps2',[filename 'a.eps'])
  figure(2)
  print('-deps2',[filename 'b.eps'])
else
  print('-deps2',[filename '.eps'])
end

%-------------------------------------------------------
%
% level l estimator
%

function sums = evsi_l(l,N)

% global option

nf = 16*2^l;
%nf = 2^l;
nc = nf/2;

    nsmp = 300;
    ns = nsmp/2;
    nn = nsmp/2;
    dresss_cost = 110;
    dressn_cost = 420;
    wtp = 30000;
    % Safety Data
    num_pat = 111; % Number of patients in observed data
    num_se = 27;   % Number of patients with side effects, given standard-of-care
    num_amb = 17;  % Number of patient with hospital care following side effect, given stardard-or-care
    num_death = 1; % Number of Deaths
    % Priors recovery
    m_1_3 =  0.45;
    v_1_3 = 0.02;
    p1_1_3 = ((1 - m_1_3)/v_1_3 - 1/ m_1_3) * (m_1_3 ^2);
    p2_1_3  = p1_1_3 * (1/m_1_3 - 1);

    m_2_3 = 0.35;
    v_2_3 = 0.02;
    p1_2_3 = ((1 - m_2_3)/v_2_3 - 1/ m_2_3) * (m_2_3 ^2);
    p2_2_3  = p1_2_3 * (1/m_2_3 - 1);

    % Probability reduction (=1-rho) of side effect given new treatment 
    % pi[2] <- rho * pi[1]
    m_rho = 0.65; % Mean 
    s_rho = 0.1; % Standard deviation
    %tau_rho = 1/s_rho^2; % Precision

    % Costs
    % Cost of ambulatory care 
    mu_amb = 2300; % Mean
    sd_amb = 90;  % Standard deviation
    m_amb = log(mu_amb) - 0.5*log(1 + sd_amb^2/mu_amb^2);
    s_amb = sqrt(log(1 + sd_amb^2/mu_amb^2));

    % Cost of hospitalization
    mu_hosp = 6500;
    sd_hosp = 980;
    m_hosp = log(mu_hosp) - 0.5*log(1 + sd_hosp^2/mu_hosp^2);
    s_hosp = sqrt(log(1 + sd_hosp^2/mu_hosp^2));

    %Cost of Death
    mu_death = 4200;
    sd_death = 560;
    m_death = log(mu_death) - 0.5*log(1 + sd_death^2/mu_death^2);
    s_death = sqrt(log(1 + sd_death^2/mu_death^2));
    %tau_death = 1/s_death ^ 2;
    % Effects
    % QALY on Chemo
    mu_chemo = 0.98;
    var_chemo = 0.001;
    p1_chemo = ((1 - mu_chemo)/var_chemo - 1/ mu_chemo) * (mu_chemo ^2);
    p2_chemo = p1_chemo * (1/mu_chemo - 1);

    % QALY on Amb
    mu_e_amb =  0.5;
    var_e_amb =  0.02;
    p1_amb = ((1 - mu_e_amb)/var_e_amb - 1/ mu_e_amb) * (mu_e_amb ^2);
    p2_amb = p1_amb * (1/mu_e_amb - 1);

    % QALY on hosp
    mu_e_hosp =  0.2;
    var_e_hosp = 0.03;
    p1_hosp = ((1 - mu_e_hosp)/var_e_hosp - 1/ mu_e_hosp) * (mu_e_hosp ^2);
    p2_hosp = p1_hosp * (1/mu_e_hosp - 1);

    % Number of patients in the population to model
    NN = 1000;

    % Time Horizon on Side Effects
    Th = 15;
    
    pi0_mu = (num_se + 1)/(num_pat + 2);
    pi0_s = (num_se + 1)*(num_pat + 1 - num_se)/((num_pat + 3)*(num_pat + 2)*(num_pat + 2));
    pi1_mu = pi0_mu*m_rho;
    pi1_s = pi0_s*0.01 + pi0_mu*pi0_mu*0.01 + 0.65*0.65*pi0_s;
    pi1_alpha = ((1 - pi1_mu)/pi1_s - 1/ pi1_mu) * (pi1_mu ^ 2);
    pi1_beta = pi1_alpha*(1 / pi1_mu - 1);
sums(1:6) = 0;

for N1 = 1:N
%   outer loop
    pi0 = betarnd(num_se + 1,num_pat - num_se + 1);
    rho = normrnd(m_rho,s_rho,1);
    pi1 = pi0*rho;
    gamma = betarnd(num_amb + 1,num_se - num_amb + 1, 1);
    gamma2 = betarnd(num_death + 1, num_se - num_amb - num_death  + 4,1);
    lambda_1_3_fix = betarnd(p1_1_3,p2_1_3,1);
    lambda_2_3_fix = betarnd(p1_2_3,p2_2_3,1);
    XAE0 = binornd(ns,pi0,1);
    XAE1 = binornd(nn,pi1,1);
    XHosp = binornd(XAE0 + XAE1,1- gamma,1);
    XDeath = binornd(XHosp,gamma2,1);
    lambda = -1/log(1-lambda_1_3_fix);
    lambda2 = -1/log(1-lambda_2_3_fix);
    XRec1 = XAE0 + XAE1 - XHosp;
    XRec2 = XHosp - XDeath;
    THC = exprnd(lambda,[1 XRec1]);
    TH = exprnd(lambda2,[1 XRec2]);
    posterior_pi0_alpha = num_se + 1 + XAE0;
    posterior_pi0_beta = num_pat - num_se + 1 + ns - XAE0;
    posterior_pi1_alpha = pi1_alpha + XAE1;
    posterior_pi1_beta = pi1_beta + nn - XAE1;
    posterior_m_1_3 = (m_1_3 + 1 - exp(-1/mean(THC)))/2;
    posterior_v_1_3 = v_1_3;
    posterior_p1_1_3 = ((1 - posterior_m_1_3)/posterior_v_1_3 - 1/ posterior_m_1_3) * (posterior_m_1_3 ^2);
    posterior_p2_1_3 = posterior_p1_1_3 * (1/posterior_m_1_3 - 1);
    posterior_m_2_3 = (m_2_3 + 1 - exp(-1/mean(TH)))/2;
    posterior_v_2_3 = v_2_3;
    posterior_p1_2_3 = ((1 - posterior_m_2_3)/posterior_v_2_3 - 1/ posterior_m_2_3) * (posterior_m_2_3 ^2);
    posterior_p2_2_3 = posterior_p1_2_3 * (1/posterior_m_2_3 - 1);
    posterior_gamma_alpha = num_amb + 1 + XAE0 + XAE1 - XHosp;
    posterior_gamma_beta = num_se - num_amb + 1 + XHosp;
    posterior_gamma2_alpha = num_death + 1 + XDeath;
    posterior_gamma2_beta = num_se - num_amb - num_death  + 4 + XHosp - XDeath;




% level l>0, with antithetic sampler
    pi_0a = betarnd(posterior_pi0_alpha,posterior_pi0_beta,nc,1);
    pi_1a = betarnd(posterior_pi1_alpha,posterior_pi1_beta,nc,1);
    in_rhoa = pi_1a ./ pi_0a;
    in_lambda_1_3_fixa = betarnd(posterior_p1_1_3,posterior_p2_1_3,nc,1);
    in_lambda_2_3_fixa = betarnd(posterior_p1_2_3,posterior_p2_2_3,nc,1);
    in_gammaa = betarnd(posterior_gamma_alpha,posterior_gamma_beta,nc,1);
    in_gamma2a = betarnd(posterior_gamma2_alpha,posterior_gamma2_beta,nc,1);
    lambda_1_2a = (1 - in_gammaa)./Th;
    lambda_1_3a = (1 - lambda_1_2a) .* in_lambda_1_3_fixa;
    lambda_1_1a = (1 - in_lambda_1_3_fixa) .* (1 - lambda_1_2a);
    lambda_2_4a = in_gamma2a ./ Th;
    lambda_2_3a = (1 - lambda_2_4a) .* in_lambda_2_3_fixa;
    lambda_2_2a = (1 - in_lambda_2_3_fixa) .* (1 - lambda_2_4a);
    
    pi_0b = betarnd(posterior_pi0_alpha,posterior_pi0_beta,nc,1);
    pi_1b = betarnd(posterior_pi1_alpha,posterior_pi1_beta,nc,1);
    in_rhob = pi_1b ./ pi_0b;
    in_lambda_1_3_fixb = betarnd(posterior_p1_1_3,posterior_p2_1_3,nc,1);
    in_lambda_2_3_fixb = betarnd(posterior_p1_2_3,posterior_p2_2_3,nc,1);
    in_gammab = betarnd(posterior_gamma_alpha,posterior_gamma_beta,nc,1);
    in_gamma2b = betarnd(posterior_gamma2_alpha,posterior_gamma2_beta,nc,1);
    lambda_1_2b = (1 - in_gammab)./Th;
    lambda_1_3b = (1 - lambda_1_2b) .* in_lambda_1_3_fixb;
    lambda_1_1b = (1 - in_lambda_1_3_fixb) .* (1 - lambda_1_2b);
    lambda_2_4b = in_gamma2b ./ Th;
    lambda_2_3b = (1 - lambda_2_4b) .* in_lambda_2_3_fixb;
    lambda_2_2b = (1 - in_lambda_2_3_fixb) .* (1 - lambda_2_4b);
    
    c_amba = lognrnd(m_amb,s_amb,nc,1);
    c_hospa = lognrnd(m_hosp,s_hosp,nc,1);
    c_deatha = lognrnd(m_death,s_death,nc,1);
    
    e_chemoa = betarnd(p1_chemo,p2_chemo,nc,1);
    e_amba = betarnd(p1_amb,p2_amb,nc,1);
    e_hospa = betarnd(p1_hosp,p2_hosp,nc,1);
    
    SE0a = binornd(NN,pi_0a,nc,1);
    SE1a = binornd(NN,pi_1a,nc,1);
    
    c_ambb = lognrnd(m_amb,s_amb,nc,1);
    c_hospb = lognrnd(m_hosp,s_hosp,nc,1);
    c_deathb = lognrnd(m_death,s_death,nc,1);
    
    e_chemob = betarnd(p1_chemo,p2_chemo,nc,1);
    e_ambb = betarnd(p1_amb,p2_amb,nc,1);
    e_hospb = betarnd(p1_hosp,p2_hosp,nc,1);
    
    SE0b = binornd(NN,pi_0b,nc,1);
    SE1b = binornd(NN,pi_1b,nc,1);
    
    q_chemo0a = (NN - SE0a) .* e_chemoa .* 16;
    q_chemo1a = (NN - SE1a) .* e_chemoa .* 16;
    MM_mata = zeros([4 4 nc]);
    LL_mata = zeros([4 4 nc]);
    LL_mata(1,1,:) = ones([nc 1]);
    LL_mata(2,2,:) = ones([nc 1]);
    LL_mata(3,3,:) = ones([nc 1]);
    LL_mata(4,4,:) = ones([nc 1]);
    MM_mata(1,1,:) = lambda_1_1a;
    MM_mata(1,2,:) = lambda_1_2a;
    MM_mata(1,3,:) = lambda_1_3a;
    MM_mata(2,2,:) = lambda_2_2a;
    MM_mata(2,3,:) = lambda_2_3a;
    MM_mata(2,4,:) = lambda_2_4a;
    MM_mata(3,3,:) = ones([nc 1]);
    MM_mata(4,4,:) = ones([nc 1]);
    q_chemo0b = (NN - SE0b) .* e_chemob .* 16;
    q_chemo1b = (NN - SE1b) .* e_chemob .* 16;
    MM_matb = zeros([4 4 nc]);
    LL_matb = zeros([4 4 nc]);
    LL_matb(1,1,:) = ones([nc 1]);
    LL_matb(2,2,:) = ones([nc 1]);
    LL_matb(3,3,:) = ones([nc 1]);
    LL_matb(4,4,:) = ones([nc 1]);
    MM_matb(1,1,:) = lambda_1_1b;
    MM_matb(1,2,:) = lambda_1_2b;
    MM_matb(1,3,:) = lambda_1_3b;
    MM_matb(2,2,:) = lambda_2_2b;
    MM_matb(2,3,:) = lambda_2_3b;
    MM_matb(2,4,:) = lambda_2_4b;
    MM_matb(3,3,:) = ones([nc 1]);
    MM_matb(4,4,:) = ones([nc 1]);
    for i = 1:nc
        for k = 1:15
                LL_mata(:,:,i) = LL_mata(:,:,i) + MM_mata(:,:,i)^k;
                LL_matb(:,:,i) = LL_matb(:,:,i) + MM_matb(:,:,i)^k;
        end
    end
    NBa = zeros(nc,2);
    NBb = zeros(nc,2);
    for i = 1:nc
           NBa(i,1) = -dresss_cost + (1/(NN*(Th + 1)))*[SE0a(i) 0 0 0]*LL_mata(:,:,i)*[wtp*e_amba(i)- c_amba(i);wtp*e_hospa(i) - c_hospa(i);wtp*e_chemoa(i);-c_deatha(i)] + wtp*q_chemo0a(i)/(NN*(Th + 1));
           NBa(i,2) = -dressn_cost + (1/(NN*(Th + 1)))*[SE1a(i) 0 0 0]*LL_mata(:,:,i)*[wtp*e_amba(i)- c_amba(i);wtp*e_hospa(i) - c_hospa(i);wtp*e_chemoa(i);-c_deatha(i)] + wtp*q_chemo1a(i)/(NN*(Th + 1));
           NBb(i,1) = -dresss_cost + (1/(NN*(Th + 1)))*[SE0b(i) 0 0 0]*LL_matb(:,:,i)*[wtp*e_ambb(i)- c_ambb(i);wtp*e_hospb(i) - c_hospb(i);wtp*e_chemob(i);-c_deathb(i)] + wtp*q_chemo0b(i)/(NN*(Th + 1));
           NBb(i,2) = -dressn_cost + (1/(NN*(Th + 1)))*[SE1b(i) 0 0 0]*LL_matb(:,:,i)*[wtp*e_ambb(i)- c_ambb(i);wtp*e_hospb(i) - c_hospb(i);wtp*e_chemob(i);-c_deathb(i)] + wtp*q_chemo1b(i)/(NN*(Th + 1)); 
    end
    pi_pera = binopdf(XAE0,ns,pi_0a);
    pi1_pera = binopdf(XAE1,nn,pi_1a);
    pi_perb = binopdf(XAE0,ns,pi_0b);
    pi1_perb = binopdf(XAE1,nn,pi_1b);
    qpi_pera = betapdf(pi_0a,posterior_pi0_alpha,posterior_pi0_beta);
    qpi_perb = betapdf(pi_0b,posterior_pi0_alpha,posterior_pi0_beta);
    qpi1_pera = betapdf(pi_1a,posterior_pi1_alpha,posterior_pi1_beta);
    qpi1_perb = betapdf(pi_1b,posterior_pi1_alpha,posterior_pi1_beta);
    qlambda_pera = betapdf(in_lambda_1_3_fixa,posterior_p1_1_3,posterior_p2_1_3);
    qlambda_perb = betapdf(in_lambda_1_3_fixb,posterior_p1_1_3,posterior_p2_1_3);
    qlambda2_pera = betapdf(in_lambda_2_3_fixa,posterior_p1_2_3,posterior_p2_2_3);
    qlambda2_perb = betapdf(in_lambda_2_3_fixb,posterior_p1_2_3,posterior_p2_2_3);
    
    pipi_pera = betapdf(pi_0a,num_se + 1,num_pat - num_se + 1);
    pipi_perb = betapdf(pi_0b,num_se + 1,num_pat - num_se + 1);
    pirho_pera = normpdf(in_rhoa,m_rho,s_rho);
    pirho_perb = normpdf(in_rhob,m_rho,s_rho);
    pilambda_pera = betapdf(in_lambda_1_3_fixa,p1_1_3,p2_1_3);
    pilambda_perb = betapdf(in_lambda_1_3_fixb,p1_1_3,p2_1_3);
    pilambda2_pera = betapdf(in_lambda_2_3_fixa,p1_2_3,p2_2_3);
    pilambda2_perb = betapdf(in_lambda_2_3_fixb,p1_2_3,p2_2_3);
    
    lambda_pera = prod(exppdf(repmat(THC,nc,1),-1./log(1 - in_lambda_1_3_fixa)),2);
    lambda_perb = prod(exppdf(repmat(THC,nc,1),-1./log(1 - in_lambda_1_3_fixb)),2);
    lambda2_pera = prod(exppdf(repmat(TH,nc,1),-1./log(1 - in_lambda_2_3_fixa)),2);
    lambda2_perb = prod(exppdf(repmat(TH,nc,1),-1./log(1 - in_lambda_2_3_fixb)),2);
    
    rhoythetaa = pi_pera.*pi1_pera.*lambda_pera.*lambda2_pera;
    rhoythetab = pi_perb.*pi1_perb.*lambda_perb.*lambda2_perb;
    
    qythetaa = qpi_pera.*qpi1_pera.*qlambda_pera.*qlambda2_pera;
    qythetab = qpi_perb.*qpi1_perb.*qlambda_perb.*qlambda2_perb;
    
    pithetaa = pipi_pera.*pirho_pera.*pilambda_pera.*pilambda2_pera;
    pithetab = pipi_perb.*pirho_perb.*pilambda_perb.*pilambda2_perb;
    
    vea = mean(rhoythetaa.*pithetaa./qythetaa);
    veb = mean(rhoythetab.*pithetab./qythetab);
    veab = (vea + veb)/2;
    
    vva = NBa .* rhoythetaa .* pithetaa ./qythetaa;
    vvb = NBb .* rhoythetab .*pithetab ./qythetab;
    
    Pf = mean((max(vva.').' + max(vvb.').'))/(2*veab) -  max(mean(vva + vvb)/(2*veab));
    dP = Pf;
    dP = dP + (1/2)*(-mean(max(vva.'))/vea - mean(max(vvb.'))/veb + max(mean(vva)/vea) + max(mean(vvb)/veb)); 
    
  

  sums(1) = sums(1) + sum(dP);
  sums(2) = sums(2) + sum(dP.^2);
  sums(3) = sums(3) + sum(dP.^3);
  sums(4) = sums(4) + sum(dP.^4);
  sums(5) = sums(5) + sum(Pf);
  sums(6) = sums(6) + sum(Pf.^2);
end