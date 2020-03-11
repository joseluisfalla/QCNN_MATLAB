function [df_dc] = FiniteDifference(d, epsilon, C_var, C_const, Permute, N, px, py, pz, batch_size, state)

C_var = permute(C_var, [3 1 2]);

C_p = C_var + epsilon * eye(length(C_var));

C_m = C_var - epsilon * eye(length(C_var));

fidelity_plus = zeros(length(C_var), 1);

fidelity_minus = zeros(length(C_var), 1);

for i = 1:length(C_var)
    
    C_plus = permute(C_p(:,i), [3 2 1]);
    
    C_minus = permute(C_m(:,i), [3 2 1]);
    
    fidelity_plus(i) = Fidelity(d, C_plus, C_const, Permute, N, px, py, pz, batch_size, state);
    
    fidelity_minus(i) = Fidelity(d, C_minus, C_const, Permute, N, px, py, pz, batch_size, state);

end

df_dc = (1 / (2 * epsilon)) * (fidelity_plus - fidelity_minus);