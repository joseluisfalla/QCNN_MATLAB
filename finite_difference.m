function [df_dc] = finite_difference(d, epsilon, C_variable, C_constant, N, px, py, pz, batch_size, state)

C_variable = permute(C_variable, [3 1 2]);

C_p = C_variable + epsilon * eye(length(C_variable));

C_m = C_variable - epsilon * eye(length(C_variable));

fidelity_plus = zeros(length(C_variable), 1);

fidelity_minus = zeros(length(C_variable), 1);

for i = 1:length(C_variable)
    
    C_plus = permute(C_p(:,i), [3 2 1]);
    
    C_minus = permute(C_m(:,i), [3 2 1]);
    
    fidelity_plus(i) = average_fidelity(d, C_plus, C_constant, N, px, py, pz, batch_size, state);
    
    fidelity_minus(i) = average_fidelity(d, C_minus, C_constant, N, px, py, pz, batch_size, state);

end

df_dc = (1 / (2 * epsilon)) * (fidelity_plus - fidelity_minus);