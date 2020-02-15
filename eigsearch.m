function xi = eigsearch(M, theta, RY, preci, theta_d_max, theta_range)
source_num = length(theta);
if 1 == source_num
    alpha_min = max(min(theta) - theta_range, -90);
    alpha_max = min(max(theta) + theta_range, 90);
    alpha_num = (alpha_max - alpha_min) / preci + 1;
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    Lambda = zeros(alpha_num * alpha_dnum, 1);%store of the correlations
    [E, D] = eig(RY);
    lambda = diag(D);
    temp = sum(lambda) * 0.95;
    accu = 0;
    for jj = M : -1 : 1
        accu = accu + lambda(jj, 1);
        if accu >= temp
            break;
        else
        end
    end
    E_noise = E(:, 1 : max(1, jj - 1));
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
            Psi = E_noise' * Psi;
            Lambda((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
%     [~, IX] = sort(Lambda);
%     Lambda(max(IX(1, 1) - 6 / preci * alpha_dnum, 1) : max(IX(1, 1) - 1, 1)) = max(Lambda);
%     Lambda(min(IX(1, 1) + 1, alpha_num * alpha_dnum) : min(IX(1, 1) + 6 / preci * alpha_dnum, alpha_num * alpha_dnum)) = max(Lambda);
    [~, IX] = sort(Lambda);
    xi = doasearch(alpha_min, alpha_dnum, preci, source_num, IX);
else
    alpha_min1 = max(min(theta) - theta_range, -90);
    alpha_max1 = min(min(theta) + theta_range, 90);
    alpha_min2 = max(max(theta) - theta_range, -90);
    alpha_max2 = min(max(theta) + theta_range, 90);
    alpha_num1 = (alpha_max1 - alpha_min1) / preci + 1;
    alpha_num2 = (alpha_max2 - alpha_min2) / preci + 1;
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    Lambda1 = zeros(alpha_num1 * alpha_dnum, 1);%store of the correlations
    Lambda2 = zeros(alpha_num2 * alpha_dnum, 1);%store of the correlations
    [E, D] = eig(RY);
    lambda = diag(D);
    temp = sum(lambda) * 0.95;
    accu = 0;
    for jj = M : -1 : 1
        accu = accu + lambda(jj, 1);
        if accu >= temp
            break;
        else
        end
    end
    E_noise = E(:, 1 : max(1, jj - 1));
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
            Psi = E_noise' * Psi;
            Lambda1((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
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
            Psi = E_noise' * Psi;
            Lambda2((ii - 1) * alpha_dnum + jj) = sum(diag(Psi' * Psi));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    xi_temp = zeros(2, 2);
    [~, IX] = sort(Lambda1);
    xi_temp(:, 1) = doasearch_2s(alpha_min1, alpha_dnum, preci, IX);
    [~, IX] = sort(Lambda2);
    xi_temp(:, 2) = doasearch_2s(alpha_min2, alpha_dnum, preci, IX);
    [~, I] = sort(xi_temp(1, :));
    xi = xi_temp(:, I);
end

