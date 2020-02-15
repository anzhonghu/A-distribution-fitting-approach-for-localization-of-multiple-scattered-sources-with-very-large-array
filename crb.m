function bound = crb(Nt, K, M, theta, theta_d, SNR)
sigma_s = 10 ^ (SNR * 0.1);%power of signal
sigma_n = 1;%power of noise
R = sigma_n * eye(M);
a = zeros(M, 1);
b = zeros(M, 1);
B = eye(M);
B_store = zeros(M * K, M);
theta = theta / 180 * pi;
theta_d = theta_d / 180 * pi;
for k = 1 : K
    for m = 1 : M
        a(m, 1) = exp(1i * (m - 1) * pi * sin(theta(k, 1)));
    end
    D = diag(a);
    for m = 1 : M
        for n = 1 : m
            B(m, n) = exp(-0.5 * (m - n) ^ 2 * (pi * cos(theta(k, 1)) * theta_d(k, 1)) ^ 2);
            B(n, m) = B(m, n);
        end
    end
    B_store((k - 1) * M + 1 : k * M, :) = B;
    R = R + sigma_s * D * B * D';
end
R_inv = inv(R);
parti_theta = zeros(M * K, M);
parti_theta_d = zeros(M * K, M);
parti_s = zeros(M * K, M);
W = zeros(M, M);
E = zeros(M, M);
for k = 1 : K
    for m = 1 : M
        b(m, 1) = (m - 1) * pi * cos(theta(k, 1));
        a(m, 1) = exp(1i * (m - 1) * pi * sin(theta(k, 1)));
    end
    Q = diag(b);
    D = diag(a);
    for m = 1 : M
        for n = 1 : m - 1
            W(m, n) = (m - n) ^ 2  * 0.5 * sin(2 * theta(k, 1)) * (theta_d(k, 1) * pi) ^ 2;
            W(n, m) = W(m, n);
            E(m, n) = -(m - n) ^ 2 * theta_d(k, 1) * (cos(theta(k, 1)) * pi) ^ 2;
            E(n, m) = E(m, n);
        end
    end
    parti_theta((k - 1) * M + 1: k * M, :) = sigma_s * (1i * Q * D * B_store((k - 1) * M + 1 : k * M, :) * D' ...
                                               - 1i * D * B_store((k - 1) * M + 1 : k * M, :) * D' * Q + D * (W .* B_store((k - 1) * M + 1 : k * M, :)) * D');
    parti_theta_d((k - 1) * M + 1: k * M, :) = sigma_s *  D * (E .* B_store((k - 1) * M + 1 : k * M, :)) * D';
    parti_s((k - 1) * M + 1: k * M, :) = D *  B_store((k - 1) * M + 1 : k * M, :) * D';
end

J_eta_eta = zeros(2 * K + 1, 2 * K + 1);
J_zita_eta = zeros(K, 2 * K + 1);
J_zita_zita = zeros(K, K);
for n = 1 : 2 * K + 1
    for m = 1 : n
        if n <= K
           J_eta_eta(n, m) =  trace(R \ parti_theta((n - 1) * M + 1: n * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :))); 
        else
            if n <= 2 * K
                if m <= K
                    J_eta_eta(n, m) =  trace(R \ parti_theta_d((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :))); 
                else
                    J_eta_eta(n, m) =  trace(R \ parti_theta_d((n - K - 1) * M + 1: (n - K) * M, :) * (R \ parti_theta_d((m - K - 1) * M + 1: (m - K) * M, :))); 
                end
            else
                if m <= K
                     J_eta_eta(n, m) = trace(R \ (R \ parti_theta((m - 1) * M + 1: m * M, :)));
                else
                    if m <= 2 * K
                        J_eta_eta(n, m) =  trace(R \ (R \ parti_theta_d((m - K - 1) * M + 1: (m - K) * M, :))); 
                    else
                        J_eta_eta(n, m) =  trace(R \ R_inv); 
                    end
                end
            end
        end
        J_eta_eta(m, n) = J_eta_eta(n, m);
    end
end
for n = 1 : K
    for m = 1 : K
        J_zita_eta(n, m) = trace(R \  parti_s((n - 1) * M + 1: n * M, :) * (R \ parti_theta((m - 1) * M + 1: m * M, :))); 
    end
    for m = K + 1 : 2 * K
        J_zita_eta(n, m) =  trace(R \  parti_s((n - 1) * M + 1: n * M, :) * (R \ parti_theta_d((m - K - 1) * M + 1: (m - K) * M, :))); 
    end
    m = 2 * K + 1;
    J_zita_eta(n, m) = trace((R \  parti_s((n - 1) * M + 1: n * M, :)) / R); 
end
for n = 1 : K
    for m = 1 : n
        J_zita_zita(n, m) = trace(R \  parti_s((n - 1) * M + 1: n * M, :) * (R \ parti_s((m - 1) * M + 1: m * M, :))); 
        J_zita_zita(m, n) = J_zita_zita(n, m);
    end
end
CRB_eta = diag(inv(J_eta_eta - (J_zita_eta.' / J_zita_zita) * J_zita_eta)) / Nt;
bound = sqrt(real(CRB_eta(1 : 2 * K))) * 180 / pi;





