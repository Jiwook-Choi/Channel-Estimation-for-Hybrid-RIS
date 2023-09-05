function [beta,gamma] = Func_Pathloss(R)
global K N N_RF L

%%% Random points in a circle
t = 2*pi*rand(K,1);
r = R*sqrt(rand(K,1));
x = r.*cos(t);
y = r.*sin(t);

%%%
d_k = sqrt(abs(x-20).^2 + abs(y-20).^2);
d_h = sqrt(abs(80).^2 + abs(20).^2);

beta = 10^(-20/10)*d_k.^-2.1;    % Pathloss for UEs-to-RIS Channel
gamma = 10^(-20/10)*d_h.^-2.2;   % Pathloss for RIS-to-BS Channel

beta = sort(beta,'descend');     % WLOG..
end