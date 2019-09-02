function [error_state_dm] = ErrorChannel(N, p_x, p_y, p_z, state)

p = p_x + p_y + p_z;

PX = [0, 1; 1, 0]; %Pauli X-Gate (bit flip)
PY = [0, -1i; 1i, 0]; %Pauli Y-Gate (bit-phase flip)
PZ = [1, 0; 0, -1]; %Pauli Z-Gate (phase flip)

X = cell(1,N);
Y = cell(1,N);
Z = cell(1,N);

state_dm = state * state';

index = 0;

for i = 1:N
    X{i} = kron(kron(eye(2^(index)), PX), eye(2^((N - 1) - index)));
    Y{i} = kron(kron(eye(2^(index)), PY), eye(2^((N - 1) - index)));
    Z{i} = kron(kron(eye(2^(index)), PZ), eye(2^((N - 1) - index)));
    
    index = index + 1;
    
end

%Z12 = kron(kron(PZ, PZ), eye(128));
%Z23 = kron(kron(kron(eye(2), PZ), PZ), eye(64));
%Z45 = kron(kron(kron(eye(8), PZ), PZ), eye(16));
%Z56 = kron(kron(kron(eye(16), PZ), PZ), eye(8));
%Z78 = kron(kron(kron(eye(64), PZ), PZ), eye(2));
%Z89 = kron(kron(eye(128), PZ), PZ);

error_state_dm = ((1 - p) * state_dm) + p_x * (X{1} * state_dm * X{1} + X2 * state_dm * X2 + X3 * state_dm * X3 + X4 * state_dm * X4 ...
    + X5 * state_dm * X5 + X6 * state_dm * X6 + X7 * state_dm * X7 + X8 * state_dm * X8 + X9 * state_dm * X9) + p_y * (Y1 * state_dm * Y1 ...
    + Y2 * state_dm * Y2 + Y3 * state_dm * Y3 + Y4 * state_dm * Y4 + Y5 * state_dm * Y5 + Y6 * state_dm * Y6 + Y7 * state_dm * Y7 ...
    + Y8 * state_dm * Y8 + Y9 * state_dm * Y9) + p_z * (Z1 * state_dm * Z1 + Z2 * state_dm * Z2 + Z3 * state_dm * Z3 + Z4 * state_dm * Z4 ...
    + Z5 * state_dm * Z5 + Z6 * state_dm * Z6 + Z7 * state_dm * Z7 + Z8 * state_dm * Z8 + Z9 * state_dm * Z9);

