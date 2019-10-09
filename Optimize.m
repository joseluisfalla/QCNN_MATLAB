function [fidelity_update, eta] = Optimize(d, epsilon, eta, C, Permute, N, px, py, pz, batch_size, iters, state)

fidelity_update = zeros(iters, 1);

for i = 1:iters

    fidelity_previous = Fidelity(d, C, Permute, N, px, py, pz, batch_size, state);

    df_dc = FiniteDifference(d, epsilon, C, Permute, N, px, py, pz, batch_size, state);

    C = permute(C, [3 1 2]);

    C_update = C + eta * df_dc;

    C_update = permute(C_update, [3 2 1]);

    fidelity_update(i) = Fidelity(d, C_update, Permute, N, px, py, pz, batch_size, state);

    if (fidelity_update - fidelity_previous) > 0

        eta = eta + (0.05 * eta);

    elseif (fidelity_update - fidelity_previous) < 0

        eta = eta - (0.5 * eta);

    end

    C = C_update;
    
end
