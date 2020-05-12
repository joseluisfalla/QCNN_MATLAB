close all
clear
clc

tic

%%
N = 9; %Number of qubits
d = 8; %Dimensions of Unitaries

px = 0.0058; %Bit-flip error probability
py = 0.002; %Bit- and phase-flip error probability
pz = 0.002; %Phase-flip error probability
% pxx = 0.0002; %Correlated X-X error probability

epsilon = 1E-4; %Finite difference constant
eta = 10; %Learning parameter

batch_size = 10000; %Number of passes through error channel

%%

%Logical qubit choices, i.e., |+/- x, y, z>
plus_x = (1 / sqrt(2)) * [1; 1];
minus_x = (1 / sqrt(2)) * [1; -1];

plus_y = (1 / sqrt(2)) * [1; 1i];
minus_y = (1 / sqrt(2)) * [1; -1i];

plus_z = [1; 0];
minus_z = [0; 1];

%Choose logical qubit
state = plus_x;

%%
%Optimization

sample_size = 100;
optimized_fidelity_1 = zeros(sample_size, 1);
optimized_fidelity_2 = zeros(sample_size, 1);
C1_optimized = zeros(2^(2 * (N / 3)) - 1, 1, sample_size);
C2_optimized = zeros(2^(2 * (N / 3)) - 1, 1, sample_size);
fidelity_per_sample = zeros(sample_size, 1, sample_size);

for i = 1:sample_size

    C1 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); %Coefficients for U1
    C2 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); %Coefficients for U2
    
    [optimized_fidelity_1(i), C1_optimized(:,:,i), fidelity_per_sample(:,:,i), eta] = optimization(d, epsilon, eta, C1, C2, N, px, py, pz, batch_size, sample_size, state);
    
end

[max_optimized_fidelity_1, max_index_1] = max(optimized_fidelity_1);
C1_optimized = permute(C1_optimized(:,:,max_index_1), [3 2 1]);

for i = 1:sample_size
    
    C2 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]);
    
    [optimized_fidelity_2(i), C2_optimized(:,:,i), fidelity_per_sample(:,:,i), eta] = optimization(d, epsilon, eta, C2, C1_optimized, N, px, py, pz, batch_size, sample_size, state);
    
end

[max_optimized_fidelity_2, max_index_2] = max(optimized_fidelity_2);
C2_optimized = permute(C2_optimized(:,:,max_index_2), [3 2 1]);

toc