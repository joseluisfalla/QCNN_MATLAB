function df_dc = finite_difference(d, epsilon, C1, N, p, batch_size, state)

C1 = permute(C1, [3 1 2]);

C1_p = C1 + epsilon * eye(length(C1));

C1_m = C1 - epsilon * eye(length(C1));

fidelity_plus = zeros(length(C1), 1);

fidelity_minus = zeros(length(C1), 1);

for i = 1:length(C1)
    
    C1_plus = permute(C1_p(:,i), [3 2 1]);
    
    C1_minus = permute(C1_m(:,1), [3 2 1]);
    
    fidelity_plus(i) = average_fidelity(d, C1_plus, N, p, batch_size, state);
    
    fidelity_minus(i) = average_fidelity(d, C1_minus, N, p, batch_size, state);
    
end

df_dc = (1 / (2 * epsilon)) * (fidelity_plus - fidelity_minus);