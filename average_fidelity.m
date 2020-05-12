function fidelity = average_fidelity(d, C1, N, p, batch_size, state)

%Generate Unitaries
GM = GellMann(d);
U1 = expm(-1i * sum(C1 .* GM, 3));

%Encoding Mechanism
encode_state = U1' * state;
encode_state_dm = encode_state * encode_state';

PX = [0, 1; 1, 0]; %Pauli-X matrix

fidelity_per_iteration = zeros(batch_size, 1);

for i = 1:batch_size

    r = rand(1);

    if p > r

        index = (randi(N, 1) - 1);
        X = kron(kron(eye(2^(index)), PX), eye(2^((N - 1) - index)));
        error_state_dm = X * encode_state_dm * X;

        decode_state_dm = U1 * error_state_dm * U1';

        fidelity_per_iteration(i) = real(encode_state' * decode_state_dm * encode_state);

    else

        fidelity_per_iteration(i) = 1;

    end

end

fidelity = mean(fidelity_per_iteration);