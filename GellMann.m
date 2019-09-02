%% Function for constructing unitaries from Generalized Gell-Mann Matrices

function [GM1, GM2] = GellMann(d)

%Dimension of Gell-Mann Matrices
GM = zeros(d, d, d^2);
k = 1;

for i = 0:(d-1)
    for j = 0:(d-1)
        GM(:, :, k) = GenGellMann(i, j, d);
        k = k + 1;
    end
end

GM1 = GM(:, :, 2:end); %Ignore Identity Matrix
GM2 = GM(:, :, 2:end);