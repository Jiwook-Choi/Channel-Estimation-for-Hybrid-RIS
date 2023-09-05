clc;
clear;
rng('shuffle');
global K N N_RF L M

N_repeat = 10; % # of experiments (iteration numbers)

K = 4;         % # of uplink UEs
N = 32;        % # of RIS elements
N_RF = 1;      % # of RF chains in RIS <= K
M = 8;         % # of BS raon Ratio, 1 - alpha : reflection ratio
L = N*K/N_RF;  % The minimum pilot length L..

Pilot_Type = 'CD';   % TD - Time division, CD - Code division
Phi_Type = 'P';      % P - proposed, R - Randomized
Psi_Type = 'P';      % P - proposed, R - Randomized

%%% Define Power and Noise Variance
dB_scale  = 0; % 30 dBm - Transmit power
noise_spectral_density = -203; % -173dBm

noise_pwr = 10^(-203/10)*10^2*10^6; % BW 100M, linear scale..
pwr = power(10,dB_scale/10);

%%% Aborption ratio range..
alpha_temp = [0.1, 0.1];
for i = 1 : 8
    alpha_temp = [alpha_temp alpha_temp(end-1)*10^-1 (alpha_temp(end-1)*10^-1)/2 ];
end
alpha_temp(1) = [];
alpha_range = [1 - flip(alpha_temp),0.8,0.7,0.6,0.5,0.4,0.3,0.2,alpha_temp ];

%%% Define vectors for MSEs
MSE_G_Sim = zeros(length(alpha_range),length(dB_scale)); % Simulated MSE for G
MSE_H_Sim = zeros(length(alpha_range),length(dB_scale)); % Simulated MSE for H

MSE_G_Ana_LS = zeros(length(alpha_range),length(dB_scale)); % Analytic LS MSE of G
MSE_H_Ana_LS = zeros(length(alpha_range),length(dB_scale)); % Analytic LS MSE of H

MSE_G_Ana_CRB = zeros(length(alpha_range),length(dB_scale)); % Analytic CRB of G
MSE_H_Ana_CRB = zeros(length(alpha_range),length(dB_scale)); % Analytic CRB of H

Normalized_term_G = zeros(length(alpha_range),length(dB_scale)); % Normalization term for NMSE of G
Normalized_term_H = zeros(length(alpha_range),length(dB_scale)); % Normalization term for NMSE of H

P_pi = Func_Permutation_Mat(); % Permutation matrix generation

for i_repeat = 1 : N_repeat
    %%% Pilot and Phase Generation
    [S,Phi,Psi] = Func_Gen_S_Phi_Psi(Pilot_Type,Phi_Type,Psi_Type);


    %%% Random Channel Generation
    [beta,gamma] = Func_Pathloss(10); % Pathloss
    G = (randn(N,K) + 1j*randn(N,K))/sqrt(2)*diag(sqrt(beta)); % RIS-UE channel

    kappa = 0;
    H_bar = sqrt(kappa/(kappa+1))*ones(N,M);
    H_tilde = sqrt(1/(kappa+1))*(randn(N,M) + 1j*randn(N,M))/sqrt(2);
    H = sqrt(gamma)*(H_bar + H_tilde); % BS-RIS channel


    w_RIS = sqrt(noise_pwr)*(randn(N,L) + 1j*randn(N,L))/sqrt(2); % noise at the RIS
    w_BS = sqrt(noise_pwr)*(randn(M,L) + 1j*randn(M,L))/sqrt(2); % noise at the BS

    for a = 1 : length(alpha_range)
        alpha = alpha_range(a); % Absoption Ratio at the RIS, 1 - alpha : reflection ratio
        for snr_ = 1 : length(dB_scale)
            [i_repeat snr_]

            %%%%% Channel Estimation Phase
            %%% Input-Ouput Relationships
            R_RIS = zeros(N,L);     % Radiation matrix of the RIS
            Y_RIS = zeros(N_RF,L);  % Received matrix of the RIS

            Y_BS = zeros(M,L);      % Received matrix of the BS

            for t = 1 : L
                R_RIS(:,t) = sqrt(pwr(snr_))*G*S(t,:).';
                Y_RIS(:,t) = sqrt(alpha/N_RF)*Phi(:,:,t)*R_RIS(:,t) + Phi(:,:,t)/sqrt(N_RF)*w_RIS(:,t); % Received signal in RIS

                x_RIS = sqrt(1-alpha)*(Psi(t,:).').* R_RIS(:,t);
                Y_BS(:,t) = H.'*x_RIS + w_BS(:,t); % Received signal at the BS
            end
            Y_RIS = Y_RIS.';
            Y_BS = Y_BS.';

            %%% G and H channel estimations
            G_hat = Func_G_Estimator(Pilot_Type,Phi_Type,Y_RIS,S,Phi,alpha,pwr(snr_),noise_pwr,P_pi);
            H_hat = Func_H_Estimator(Pilot_Type,Psi_Type,Y_BS,G_hat,S,Psi,alpha,pwr(snr_));

            %%% Sim. MSE Calculations
            MSE_G_Sim(a,snr_) = MSE_G_Sim(a,snr_) + sum(sum(abs(G_hat - G).^2)); % Normalized MSE
            MSE_H_Sim(a,snr_) = MSE_H_Sim(a,snr_) + sum(sum(abs(H_hat - H).^2));

            Normalized_term_G(a,snr_) = Normalized_term_G(a,snr_) + sum(sum(abs(G).^2));
            Normalized_term_H(a,snr_) = Normalized_term_H(a,snr_) + sum(sum(abs(H).^2));

            %%% Analytic Results
%             MSE_G_Ana_CRB(a) = MSE_G_Ana_CRB(a) + (N*N*K)/(L*N_RF*alpha*pwr(snr_))*noise_pwr;
%             MSE_G_Ana_LS(a) = MSE_G_Ana_LS(a) + (N)/(alpha*pwr(snr_))*noise_pwr;
% 
%             temp_CRB = 0;
%             temp_LS = 0;
%             for k = 1: K
%                 beta_temp = beta;
%                 beta_temp(k) = 0;
%                 beta_temp = ones(K,1) - beta_temp./(beta(k));
%                 temp_CRB = temp_CRB + (-1)^(k-1)*log(beta(k))/(beta(k)*prod(abs(beta_temp)))*noise_pwr;
%             end
%             %%% Ana. MSE Calculations
%             MSE_H_Ana_CRB(a) = MSE_H_Ana_CRB(a) + M*N/(L*(1-alpha))*temp_CRB;
%             MSE_H_Ana_LS(a) = MSE_H_Ana_LS(a) + M*N/(pwr(snr_)*K)*(gamma/alpha + 1/(N*(1-alpha)))*temp_CRB;

        end
    end
end

Normalized_term_G = Normalized_term_G/N_repeat;
Normalized_term_H = Normalized_term_H/N_repeat;

NMSE_G = (MSE_G_Sim/N_repeat)./Normalized_term_G;
NMSE_H = (MSE_H_Sim/N_repeat)./Normalized_term_H;

NMSE_G_Ana_CRB = (MSE_G_Ana_CRB/N_repeat)./Normalized_term_G;
NMSE_H_Ana_CRB = (MSE_H_Ana_CRB/N_repeat)./Normalized_term_H;

NMSE_G_Ana_LS = (MSE_G_Ana_LS/N_repeat)./Normalized_term_G;
NMSE_H_Ana_LS = (MSE_H_Ana_LS/N_repeat)./Normalized_term_H;

loglog(NMSE_G,NMSE_H,'-x');  hold on;

% loglog(NMSE_G_Ana_CRB,NMSE_H_Ana_CRB,'-x');  hold on; % For L = NK
% loglog(NMSE_G,NMSE_H,'-x');  hold on;   % For L = NK

xlabel('NMSE of G');
xlabel('NMSE of H');

axis([10^-10 10^0 10^-5 10^0]) %[xmin xmax ymin ymax]

