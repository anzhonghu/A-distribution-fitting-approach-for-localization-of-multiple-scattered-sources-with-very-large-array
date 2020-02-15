function gamma = comet(M, theta, RY, preci, theta_range)
alpha_min = max(min(theta) - theta_range, -90);
alpha_max = min(max(theta) + theta_range, 90);
alpha_num = (alpha_max - alpha_min) / preci + 1;
Lambda = zeros(alpha_num, 1);%store
a = zeros(M, 1);
H = zeros(M, 1);
J = zeros(M * M, M);
for n = 1 : M
    for l = 1 : M
        for k = 1 : M
            if 0 == abs(l - n) - (k - 1)
                J((n - 1) * M + l, k) = 1;
            else
            end
        end
    end
end
Y = zeros(M, M);
Y(1, 1) = M;
for m = 1 : M - 1
    Y(m + 1, m + 1) = 2 * (M - m);
end
ry = reshape(RY, M * M, 1);
for ii = 1 : alpha_num
    alpha = (alpha_min + (ii - 1) * preci) / 180 * pi;
    for m = 1 : M
        a(m, 1) = exp(1i * pi * (m - 1) * sin(alpha));
    end
    Phi = diag(a);
    Psi = kron(Phi', Phi);
    y = J.' * Psi' * ry;
    beta = real(Y \ y);
    for m = 1 : M
        if beta(m, 1) >= 0
            H(m, 1) = 1;
        else
            H(m, 1) = 0;
        end
    end
    MU = diag(H ./ diag(Y));
    Lambda(ii) = ry' * ry - y.' * MU * y;
end
%%%%%%%%%%%%%%%%%%%%%%%%%
[~, I] = min(real(Lambda));
gamma = alpha_min + (I - 1) * preci;

