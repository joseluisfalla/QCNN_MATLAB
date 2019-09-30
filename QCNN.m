close all
clear
clc
tic
%%
N = 9; % Number of qubits
d = 8; % Dimensions of Gell-Mann matrices

px = 0.006;
py = 0.002;
pz = 0.002;
%pxx = 0.002;

epsilon = 1E-4; % Step for finite-difference method
eta = 10; % Learning rate

batch_size = 1000;
iters = 100; % Number of iterations for gradient descent

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

fidelity_current = Optimize(d, epsilon, eta, C1, Permute, N, px, py, pz, batch_size, iters, state);

toc
