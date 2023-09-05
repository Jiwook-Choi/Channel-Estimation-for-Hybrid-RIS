function [H_hat] = Func_H_Estimator(Pilot_Type,Psi_Type,Y_BS,G_hat,S,Psi,alpha,pwr)
global K N N_RF L M

switch Pilot_Type
    case 'TD'
        X_RIS =  Psi.*(S*G_hat.');
        H_hat = pinv(X_RIS)*Y_BS/sqrt((1-alpha)*pwr);
    case 'CD'
        switch Psi_Type
            case 'P'
                if(N_RF == 1)
                    X_RIS =  Psi.*(S*G_hat.');
                    D_hat = diag(sum(abs(G_hat).^2,2));
                    H_hat = 1/(N*K*sqrt((1-alpha)*pwr))*inv(D_hat)*X_RIS'*Y_BS;
                else
                    X_RIS =  sqrt((1-alpha)*pwr)*Psi.*(S*G_hat.');
                    H_hat = pinv(X_RIS)*Y_BS;
                end
            case 'R'
                X_RIS =  Psi.*(S*G_hat.');
                H_hat = pinv(X_RIS)*Y_BS/sqrt((1-alpha)*pwr);
        end
end
end