close all
clear
clc
tic
%%
N = 9; % Number of qubits
d = 8; % Dimensions of Gell-Mann matrices

px = 0.058;
py = 0.02;
pz = 0.02;
pxx = 0.002;

epsilon = 1E-4; % Step for finite-difference method
eta = 10; % Learning rate

%%

C1 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); % Coefficients for U1
C2 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); % Coefficients for U2

%%

SS = kron(SwapOperator(2), eye(2)) * kron(eye(2), SwapOperator(2)) * kron(SwapOperator(2), eye(2));
Permute = kron(kron(kron(kron(eye(2), SS), eye(2)), SS), eye(2));

%%

plus_x = (1 / sqrt(2)) * [1; 1];
minus_x = (1 / sqrt(2)) * [1; -1];

plus_y = (1 / sqrt(2)) * [1; 1i];
minus_y = (1 / sqrt(2)) * [1; -1i];

plus_z = [1; 0];
minus_z = [0; 1];

state = plus_x; % Change as necessary

%%

iters = 100;

fid_per_iter = zeros(1, iters);

for i = 1:iters

    fidelity_initial = Fidelity(d, C1, Permute, N, px, py, pz, pxx, state);

    fidelity_plus = Fidelity(d, C1 + epsilon, Permute, N, px, py, pz, pxx, state);

    fidelity_minus = Fidelity(d, C1 - epsilon, Permute, N, px, py, pz, pxx, state);

    dfdc = (1 / (2 * epsilon)) * (fidelity_plus - fidelity_minus);

    C1_update = C1 + eta * dfdc;

    fidelity_current = Fidelity(d, C1_update, Permute, N, px, py, pz, pxx, state);

    if (fidelity_current - fidelity_initial) > 0
        eta = eta - 0.05 * eta;
        C1 = C1_update;
    elseif (fidelity_current - fidelity_initial) < 0
        eta = eta + 0.5 * eta;
        C1_update = 0;
    end
    
    fid_per_iter(i) = fidelity_current;
    
end

toc
