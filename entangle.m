function [state] = entangle(ancilla, state_ket)

state = kron(kron(kron(kron(kron(kron(kron(kron(ancilla, ancilla), ancilla), ancilla), state_ket), ancilla), ancilla), ancilla), ancilla);