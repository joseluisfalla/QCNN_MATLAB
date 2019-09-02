function [encode_state] = Encode(d, C, Permute, state)

GM1 = GellMann(d);

U1 = expm(-1i * sum(C .* GM1, 3));
U2 = eye(d);

state = Entangle(state);

U2_encode = kron(kron(eye(d), U2'), eye(d));
U1_encode = kron(kron(U1', U1'), U1');


encode_state = U1_encode * ((Permute * (U2_encode * state)));