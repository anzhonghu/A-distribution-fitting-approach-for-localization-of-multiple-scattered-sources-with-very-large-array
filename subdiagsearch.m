function psi = subdiagsearch(M, theta, RY, preci, theta_range)
alpha_min = max(min(theta) - theta_range, -90);
alpha_max = min(max(theta) + theta_range, 90);
alpha_num = (alpha_max - alpha_min) / preci + 1;
Lambda = zeros(alpha_num, 1);%store of the values
z = zeros(M - 1, 1);
for m = 1 : M - 1
    for n = 1 : M - m
        z(m, 1) = z(m, 1) + RY(m + n, n);
    end
end
for ii = 1 : alpha_num
    alpha = (alpha_min + (ii - 1) * preci) / 180 * pi;
    for k = 1 : M - 1
        Lambda(ii, 1) = Lambda(ii, 1) + (M - k) * conj(z(k, 1)) * exp(1i * k * pi * sin(alpha));
    end
end
[~, I] = max(real(Lambda));
psi = alpha_min + (I - 1) * preci;

