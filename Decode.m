function [decode_state] = Decode(d, C1, C2, Permute, state_dm)

GM = GellMann(d);

U1 = expm(-1i * sum(C1 .* GM, 3));
U2 = expm(-1i * sum(C2 .* GM, 3));

U1_decode = kron(kron(U1, U1), U1);
U2_decode = kron(kron(eye(d), U2), eye(d));

decode_state = U2_decode * Permute * (U1_decode * state_dm * U1_decode') * Permute' * U2_decode';