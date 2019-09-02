function [fid] = Fidelity(d, C1, Permute, N, px, py, pz, pxx, input_state)

num_iters = 1000;

fpi = zeros(num_iters, 1);

entangle_state = Entangle(input_state);

for i = 1:num_iters

    error_state_dm = StochasticErrorChannel(d, C1, Permute, N, px, py, pz, pxx, input_state);

    decode_state_dm = Decode(d, C1, Permute, error_state_dm);

    fpi(i) = real(entangle_state' * decode_state_dm * entangle_state);
    
end

fid = sum(fpi) / num_iters;