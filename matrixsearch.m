function phi = matrixsearch(M, theta, RY, preci, theta_d_max, theta_range)
source_num = length(theta);
if 1 == source_num
    alpha_min = max(min(theta) - theta_range, -90);
    alpha_max = min(max(theta) + theta_range, 90);
    alpha_num = (alpha_max - alpha_min) / preci + 1;
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    Lambda_gc = zeros(alpha_num * alpha_dnum, 1);%store of the eigenvalues
    Lambda_sub = zeros(alpha_num * alpha_dnum, 1);%store of the norms
    
    for ii = 1 : alpha_num
        alpha = alpha_min + (ii - 1) * preci;
        for jj = 1 : alpha_dnum
            alpha_d = jj * preci;
            Psi = eye(M);
            for k = 1 : M
                for l = 1 : k - 1
                    Psi(k, l) = exp(1i * pi * (k - l) * sin(alpha / 180 * pi)) * exp(-0.5 * (pi * (k - l) * (alpha_d / 180 * pi) * cos(alpha / 180 * pi))^2);%gaussian
                    Psi(l, k) = conj(Psi(k, l));
                end
            end
            Psi = RY \ Psi;
            Lambda_gc((ii - 1) * alpha_dnum + jj) = max(real(eig(Psi)));
            Lambda_sub((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
        end
    end
    %GC approach
    [~, IX] = sort(Lambda_gc);
    phi(:, 1) = doasearch(alpha_min, alpha_dnum, preci, source_num, IX);
    %Subspace approach
    [~, IX] = sort(Lambda_sub);
    phi(:, 2) = doasearch(alpha_min, alpha_dnum, preci, source_num, IX);
else
    alpha_min1 = max(min(theta) - theta_range, -90);
    alpha_max1 = min(min(theta) + theta_range, 90);
    alpha_min2 = max(max(theta) - theta_range, -90);
    alpha_max2 = min(max(theta) + theta_range, 90);
    alpha_num1 = (alpha_max1 - alpha_min1) / preci + 1;
    alpha_num2 = (alpha_max2 - alpha_min2) / preci + 1;
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    Lambda_gc1 = zeros(alpha_num1 * alpha_dnum, 1);%store of the eigenvalues
    Lambda_gc2 = zeros(alpha_num2 * alpha_dnum, 1);%store of the eigenvalues
    Lambda_sub1 = zeros(alpha_num1 * alpha_dnum, 1);%store of the norms
    Lambda_sub2 = zeros(alpha_num2 * alpha_dnum, 1);%store of the norms
    for ii = 1 : alpha_num1
        alpha = alpha_min1 + (ii - 1) * preci;
        for jj = 1 : alpha_dnum
            alpha_d = jj * preci;
            Psi = eye(M);
            for k = 1 : M
                for l = 1 : k - 1
                    Psi(k, l) = exp(1i * pi * (k - l) * sin(alpha / 180 * pi)) * exp(-0.5 * (pi * (k - l) * (alpha_d / 180 * pi) * cos(alpha / 180 * pi))^2);%gaussian
                    Psi(l, k) = conj(Psi(k, l));
                end
            end
            Psi = RY \ Psi;
            Lambda_gc1((ii - 1) * alpha_dnum + jj) = max(real(eig(Psi)));
            Lambda_sub1((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
        end
    end
    for ii = 1 : alpha_num2
        alpha = alpha_min2 + (ii - 1) * preci;
        for jj = 1 : alpha_dnum
            alpha_d = jj * preci;
            Psi = eye(M);
            for k = 1 : M
                for l = 1 : k - 1
                    Psi(k, l) = exp(1i * pi * (k - l) * sin(alpha / 180 * pi)) * exp(-0.5 * (pi * (k - l) * (alpha_d / 180 * pi) * cos(alpha / 180 * pi))^2);%gaussian
                    Psi(l, k) = conj(Psi(k, l));
                end
            end
            Psi = RY \ Psi;
            Lambda_gc2((ii - 1) * alpha_dnum + jj) = max(real(eig(Psi)));
            Lambda_sub2((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
        end
    end
    %GC approach
    phi_temp = zeros(2, 2);
    [~, IX] = sort(Lambda_gc1);
    phi_temp(:, 1) = doasearch_2s(alpha_min1, alpha_dnum, preci, IX);
    [~, IX] = sort(Lambda_gc2);
    phi_temp(:, 2) = doasearch_2s(alpha_min2, alpha_dnum, preci, IX);
    [~, I] = sort(phi_temp(1, :));
    phi(:, 1 : 2) = phi_temp(:, I);
    %Subspace approach
    [~, IX] = sort(Lambda_sub1);
    phi_temp(:, 1) = doasearch_2s(alpha_min1, alpha_dnum, preci, IX);
    [~, IX] = sort(Lambda_sub2);
    phi_temp(:, 2) = doasearch_2s(alpha_min2, alpha_dnum, preci, IX);
    [~, I] = sort(phi_temp(1, :));
    phi(:, 3 : 4) = phi_temp(:, I);
end

