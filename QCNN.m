close all
clear
clc
tic
%%
N = 9; % Number of qubits
d = 8; % Dimensions of Gell-Mann matrices

px = 0.0058;
py = 0.002;
pz = 0.002;
pxx = 0.0002;

epsilon = 1E-4; % Step for finite-difference method
eta = 10; % Learning rate

batch_size = 10000; % Number of passes through error channel

%%

SkipSwap = kron(SwapOperator(2), eye(2)) * kron(eye(2), SwapOperator(2)) * kron(SwapOperator(2), eye(2));
Permute = kron(kron(kron(kron(eye(2), SkipSwap), eye(2)), SkipSwap), eye(2));

%%

plus_x = (1 / sqrt(2)) * [1; 1];
minus_x = (1 / sqrt(2)) * [1; -1];

plus_y = (1 / sqrt(2)) * [1; 1i];
minus_y = (1 / sqrt(2)) * [1; -1i];

plus_z = [1; 0];
minus_z = [0; 1];

state = plus_x; % Change as necessary

%%

sample_size = 100;
optimized_fidelity = zeros(sample_size, 1);
C_optimized = zeros(2^(2 * (N / 3)) - 1, 1, sample_size);
fidelity_per_sample = zeros(sample_size, 1, sample_size);

for i = 1:sample_size

    C1 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); % Coefficients for U1
    C2 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); % Coefficients for U2
    
    [optimized_fidelity(i), C_optimized(:,:,i), fidelity_per_sample(:,:,i), eta] = Optimize(d, epsilon, eta, C1, C2, Permute, N, px, py, pz, batch_size, state);
    
end

[max_optimized_fidelity, max_index] = max(optimized_fidelity);
C_optimized = permute(C_optimized(:,:,max_index), [3 2 1]);

toc