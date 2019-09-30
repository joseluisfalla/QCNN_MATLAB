function [fidelity] = Fidelity(d, C, Permute, N, px, py, pz, batch_size, state)

encode_state = Encode(d, C, Permute, state);
encode_state_dm = encode_state * encode_state';

p = px + py + pz;

PX = [0, 1; 1, 0];
PY = [0, -1i; 1i, 0];
PZ = [1, 0; 0, -1];

X = cell(1, N);
Y = cell(1, N);
Z = cell(1, N);

fidelity_per_iteration = zeros(batch_size, 1);

for i = 1:batch_size
    
    r = rand(1);

    if p > r
        rx = rand(1);
        ry = rand(1);
        rz = rand(1);

        index = 0;

        for j = 1:N
            X{j} = kron(kron(eye(2^(index)), PX), eye(2^((N - 1) - index)));
            Y{j} = kron(kron(eye(2^(index)), PY), eye(2^((N - 1) - index)));
            Z{j} = kron(kron(eye(2^(index)), PZ), eye(2^((N - 1) - index)));

            if px > rx
                error_state_dm = X{j} * encode_state_dm * X{j};
                encode_state_dm = error_state_dm;
            elseif py > ry
                error_state_dm = Y{j} * encode_state_dm * Y{j};
                encode_state_dm = error_state_dm;
            elseif pz > rz
                error_state_dm = Z{j} * encode_state_dm * Z{j};
                encode_state_dm = error_state_dm;
            end

            index = index + 1;
            error_state_dm = encode_state_dm;

        end

        decode_state_dm = Decode(d, C, Permute, error_state_dm);
        fidelity_per_iteration(i) = real(encode_state' * decode_state_dm * encode_state);

    else

        fidelity_per_iteration(i) = 1;

    end
 
end

fidelity = sum(fidelity_per_iteration) / batch_size;
