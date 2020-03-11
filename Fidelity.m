function [fidelity] = Fidelity(d, C1, C2, Permute, N, px, py, pz, batch_size, state)

encode_state = Encode(d, C1, C2, Permute, state);
encode_state_dm = encode_state * encode_state';

p = px + py + pz;

PX = [0, 1; 1, 0];
PY = [0, -1i; 1i, 0];
PZ = [1, 0; 0, -1];

fidelity_per_iteration = zeros(batch_size, 1);

for i = 1:batch_size
    
    r = rand(1);

    if p > r
        
        rx = rand(1);
        ry = rand(1);
        rz = rand(1);

        if px > rx

            index = (randi(N, 1) - 1);
            X = kron(kron(eye(2^(index)), PX), eye(2^((N - 1) - index)));
            error_state_dm = X * encode_state_dm * X;

        elseif py > ry

            index = (randi(N, 1) - 1);
            Y = kron(kron(eye(2^(index)), PY), eye(2^((N - 1) - index)));
            error_state_dm = Y * encode_state_dm * Y;

        elseif pz > rz

            index = (randi(N, 1) - 1);
            Z = kron(kron(eye(2^(index)), PZ), eye(2^((N - 1) - index)));
            error_state_dm = Z * encode_state_dm * Z;
            
        else
            
            error_state_dm = encode_state_dm;

        end

        decode_state_dm = Decode(d, C1, C2, Permute, error_state_dm);
        fidelity_per_iteration(i) = real(encode_state' * decode_state_dm * encode_state);

    else

        fidelity_per_iteration(i) = 1;

    end
 
end

fidelity = sum(fidelity_per_iteration) / batch_size;