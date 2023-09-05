function [S,Phi,Psi] = Func_Gen_S_Phi_Psi(Pilot_Type,Phi_Type,Psi_Type)
global K N N_RF L

S = zeros(L,K);          % Pilot Matrix
Phi = zeros(N_RF, N, L); % Phase-shift matrix for Sensing
Psi = zeros(L,N);        % Phase-shift matrix for Reflection

F_NK = fft(eye(N*K)); % NK by NK DFT

switch Pilot_Type
    case 'TD' % Time - Division
        L_K = L/K;    % Each UE L length
        F_N = fft(eye(N));

        %%% Pilot Matrix
        for i = 1 : K
            temp_vec = zeros(1,L);
            temp_vec((i-1)*L_K+1: i*L_K) = sqrt(K);
            S(:,i) = temp_vec;
        end
        %
        %%% Analog combiner for RF in RIS
        switch Phi_Type
            case 'P'
                for r = 1 : N_RF
                    for n = 1 : N
                        for k = 1: K
                            I_k = (k-1)*L_K+1: k*L_K;
                            I_r = (r-1)*L_K+1: r*L_K;
                            Phi(r,n,I_k) = F_N(I_r,n);
                        end
                    end
                end
            case 'R'
                for r = 1 : N_RF
                    for n = 1 : N
                        Phi(r,n,:) = exp(2*pi*1j*rand(1,L));
                    end
                end
        end
        % %%% Reflection phase at the RIS
        F_Lk = fft(eye(L_K));
        for n = 1 : N
            switch Psi_Type
                case 'P'
                    for k = 1 : K
                        I_k = (k-1)*L_K+1: k*L_K;
                        Psi(I_k,n) = F_Lk(:, mod(n-1,L_K)+1);
                    end
                case 'R'
                    Psi(:,n) = exp(2*pi*1j*rand(1,L));
            end
        end

    case 'CD' % Code - Division
        %%% Pilot Matrix
        for k = 1 : K
            S(:,k) = F_NK(1:L,(k-1)*N+1); % N interval -> Orthogonal sequences..
        end


        %%% Analog combiner for RF in RIS
        switch Phi_Type
            case 'P'
                for r = 1 : N_RF
                    I_r = (r-1)*L + 1: r*L;
                    for n= 1 : N
                        Phi(r,n,:) = F_NK(I_r,n);
                    end
                end
            case 'R'
                Phi = exp(2*pi*1j*rand(N_RF,N,L));
        end


        %%% Reflection phase at the RIS
        for n = 1 : N
            switch Psi_Type
                case 'P'
                    Psi(:,n) = F_NK(1: L, (n-1)*N_RF + 1);
                case 'R'
                    Psi(:,n) = exp(2*pi*1j*rand(1,L));
            end
        end
end
end