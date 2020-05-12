function [state] = entangle(state_ket)

state_zero = [1; 0];

state = kron(kron(kron(kron(kron(kron(kron(kron(state_zero, state_zero), state_zero), state_zero), state_ket), state_zero), state_zero), state_zero), state_zero);