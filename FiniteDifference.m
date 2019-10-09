function [df_dc] = FiniteDifference(d, epsilon, C, Permute, N, px, py, pz, batch_size, state)

C = permute(C, [3 1 2]);

C_p = C + epsilon * eye(length(C));

C_m = C - epsilon * eye(length(C));

fidelity_plus = zeros(length(C), 1);

fidelity_minus = zeros(length(C), 1);

for i = 1:length(C)
    
    C_plus = permute(C_p(:,i), [3 2 1]);
    
    C_minus = permute(C_m(:,i), [3 2 1]);
    
    fidelity_plus(i) = Fidelity(d, C_plus, Permute, N, px, py, pz, batch_size, state);
    
    fidelity_minus(i) = Fidelity(d, C_minus, Permute, N, px, py, pz, batch_size, state);

end

df_dc = (1 / (2 * epsilon)) * (fidelity_plus - fidelity_minus);
