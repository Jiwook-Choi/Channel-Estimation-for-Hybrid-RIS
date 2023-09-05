function [P_pi] = Func_Permutation_Mat()
global K N N_RF L

vec_temp = 1:L*N_RF;

Y = reshape(vec_temp, L, N_RF);

Y_vec = reshape(Y,[],1);
Y_t_vec = reshape(Y',[],1);

P_pi = zeros(L*N_RF, L*N_RF);

for i = 1 : L*N_RF
    idx = find(Y_t_vec == Y_vec(i));
    P_pi(i,idx) = 1;
end
end
