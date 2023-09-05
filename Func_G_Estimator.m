function [G_hat] = Func_G_Estimator(Pilot_Type,Phi_Type,Y_RIS,S,Phi,alpha,pwr,noise_pwr,P_pi)
global K N N_RF L

switch Pilot_Type
    %%%% Time-Division
    case 'TD'
        L_K = L/K;    % Each UE L length
        F_N = fft(eye(N)); % N by N DFT

        switch Phi_Type
            case 'P'
                G_hat = zeros(N,K);
                for k = 1 : K
                    I_k = (k-1)*L_K+1: k*L_K;
                    Omega_k = zeros(N,N);
                    for n = 1 : N
                        Omega_k(:,n) = sqrt(K)*reshape(squeeze(Phi(:,n,I_k)).',N,1);
                    end
                    vec_y_RIS_k = reshape(Y_RIS(I_k,:),[],1);
                    G_hat(:,k) = 1/N*sqrt(N_RF/(alpha*pwr*K))*F_N'*vec_y_RIS_k;
                end
            case 'R'
                vec_y_RIS = reshape(Y_RIS,[],1); %vectorization
                %%% Omeaga matrix construction
                bar_Phi = [];
                for n = 1 : N_RF
                    bar_Phi = [bar_Phi; squeeze(Phi(n,:,:)).'  ];
                end
                Omega = kron(ones(1,K),bar_Phi) .* kron(kron(ones(N_RF,1), S),ones(1,N));

                %%%  Covariance matrix of W_RIS
                temp = [];
                for t = 1 : L
                    temp = blkdiag(temp,noise_pwr/N_RF*Phi(:,:,t)*Phi(:,:,t)');
                end
                C = P_pi*temp*P_pi.';
                C_inv = inv(C);

                %%% Generalized LS..
                g_hat = inv(Omega'*C_inv*Omega)*Omega'*C_inv*vec_y_RIS*sqrt(N_RF/(alpha*pwr));
                G_hat = reshape(g_hat,N,K);
        end

        %%%% Code-Division
    case 'CD'
        vec_y_RIS = reshape(Y_RIS,[],1);
        switch Phi_Type
            case 'P'
                F_NK = fft(eye(N*K)); % NK by NK DFT
                g_hat = (F_NK'*vec_y_RIS)*1/(N*K)*sqrt(N_RF/(alpha*pwr));

            case 'R'
                %%% Omeaga matrix construction
                bar_Phi = [];
                for n = 1 : N_RF
                    bar_Phi = [bar_Phi; squeeze(Phi(n,:,:)).'  ];
                end
                Omega = kron(ones(1,K),bar_Phi) .* kron(kron(ones(N_RF,1), S),ones(1,N));

                %%%  Covariance matrix of W_RIS
                temp = [];
                for t = 1 : L
                    temp = blkdiag(temp,noise_pwr/N_RF*Phi(:,:,t)*Phi(:,:,t)');
                end
                C = P_pi*temp*P_pi.';
                C_inv = inv(C);

                %%% Generalized LS..
                g_hat = inv(Omega'*C_inv*Omega)*Omega'*C_inv*vec_y_RIS*sqrt(N_RF/(alpha*pwr));
        end
        G_hat = reshape(g_hat,N,K);
end
end