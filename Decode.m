function [decode_state] = Decode(d, C, Permute, state_dm)

GM1 = GellMann(d);

U1 = expm(-1i * sum(C .* GM1, 3));
U2 = eye(d);

U1_decode = kron(kron(U1, U1), U1);
U2_decode = kron(kron(eye(d), U2), eye(d));

decode_state = U2_decode * Permute * (U1_decode * state_dm * U1_decode') * Permute' * U2_decode';