close all
clear
clc

tic

%%

N = 3; %Set number of qubits
d = 8; %Set dimension of Unitaries

p = 0.01; %Set error probability (bit-flip)

epsilon = 1E-4; %Finite-difference parameter
eta = 10; %Learning parameter

batch_size = 100000; %Number of runs to compute average fidelity

zero_ket = [1; 0]; %Ancilla state ket |0>
plus_x_ket = (1 / sqrt(2)) * [1; 1]; %Logical qubit to be encoded |+x>

state = kron(kron(zero_ket, plus_x_ket), zero_ket); %Product state

sample_size = 100; %Number of optimization iterations

optimized_fidelity = zeros(sample_size, 1); %Initialize optimized fidelity vector
C_optimized = zeros(2^(2 * N) - 1, 1 , sample_size); %Initialize optimized coefficient matrix
fidelity_per_sample = zeros(sample_size, 1, sample_size);

%%
%Optimization

for i = 1:sample_size
    
    C1 = permute(((2 * pi) .* rand((d^2 - 1), 1)), [3 2 1]); %Initial random coefficient vector {0, 2pi}

    [optimized_fidelity(i), C_optimized(:,:,i), fidelity_per_sample(:,:,i), eta] = optimization(d, epsilon, eta, C1, N, p, batch_size, sample_size, state);
    
end

[max_optimized_fidelity, max_index] = max(optimized_fidelity); %Find maximum value of fidelity
C_optimized = permute(C_optimized(:,:,max_index), [3 2 1]); %Find corresponding coefficient vector

%%
%Test optimized coefficients against a random vector

fid_rand = zeros(sample_size, 1);
fid_opt = zeros(sample_size, 1);

for j = 1:sample_size
    
    fid_rand(j,1) = average_fidelity(d, C1, N, p, batch_size, state);
    
    fid_opt(j,1) = average_fidelity(d, C_optimized, N, p, batch_size, state);
    
end

plot(1:sample_size, fid_rand, 'r+:', 1:sample_size, fid_opt, 'b*:');

toc