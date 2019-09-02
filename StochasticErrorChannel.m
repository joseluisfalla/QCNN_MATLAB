function [error_state_dm] = StochasticErrorChannel(d, C1, Permute, N, px, py, pz, pxx, state)

state = Encode(d, C1, Permute, state);

state_dm = state * state';

p = px + py + pz;
r = rand(1);
QS = [1,2,4,5,7,8]; % Qubit indeces for correlated error model.

PX = [0, 1; 1, 0];
PY = [0, -1i; 1i, 0];
PZ = [1, 0; 0, -1];

X = cell(1,N);
Y = cell(1,N);
Z = cell(1,N);

index = 0;

for i = 1:N
    X{i} = kron(kron(eye(2^(index)), PX), eye(2^((N - 1) - index)));
    Y{i} = kron(kron(eye(2^(index)), PY), eye(2^((N - 1) - index)));
    Z{i} = kron(kron(eye(2^(index)), PZ), eye(2^((N - 1) - index)));

    if r < p
        rx = rand(1);
        ry = rand(1);
        rz = rand(1);
        if rx < px
            error_state_dm = X{i} * state_dm * X{i};
            state_dm = error_state_dm;
        elseif ry < py
            error_state_dm = Y{i} * state_dm * Y{i};
            state_dm = error_state_dm;
        elseif rz < pz
            error_state_dm = Z{i} * state_dm * Z{i};
            state_dm = error_state_dm;
        end
    end
    
    error_state_dm = state_dm;
    index = index + 1;

end

for i = QS
    rxx = rand(1);
    if rxx < pxx
        error_state_dm = X{i} * X{i+1} * error_state_dm * X{i} * X{i+1};
        state_dm = error_state_dm;
    end
    
    error_state_dm = state_dm;
    
end