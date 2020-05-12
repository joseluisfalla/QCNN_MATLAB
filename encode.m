function [encode_state] = encode(d, C1, C2, Permute, state)

GM = GellMann(d);

U1 = expm(-1i * sum(C1 .* GM, 3));
U2 = expm(-1i * sum(C2 .* GM, 3));

state = entangle(state);

U2_encode = kron(kron(eye(d), U2'), eye(d));
U1_encode = kron(kron(U1', U1'), U1');


encode_state = U1_encode * ((Permute * (U2_encode * state)));