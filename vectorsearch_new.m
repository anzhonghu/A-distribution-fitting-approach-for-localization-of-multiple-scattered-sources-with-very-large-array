function eta = vectorsearch_new(Nt, M, theta, Y, preci, theta_d_max, theta_range)
source_num = length(theta);
if 1 == source_num
%     Eta = zeros(2, source_num);
%     z = zeros(M, 1);
%     a = zeros(M, 1);
%     b = zeros(M, 1);
%     alpha_min = max(min(theta) - theta_range, -90);
%     alpha_max = min(max(theta) + theta_range, 90);
%     alpha_num = (alpha_max - alpha_min) / preci + 1;
%     alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
%     len = alpha_num * alpha_dnum;
%     Lambda = zeros(len, 1);%store of the eigenvalues
%     B = zeros(M, len);
%     for ii = 1 : alpha_num
%         alpha = (alpha_min + (ii - 1) * preci) / 180 * pi;
%         for jj = 1 : alpha_dnum
%             alpha_d = jj * preci / 180 * pi;
%             for m = 1 : M
%                 theta_m = -pi * 0.5 + pi * m / M;
%                 b(m, 1) = 1 / sqrt(2 * pi) / alpha_d * exp(-0.5 * (theta_m - alpha) ^ 2 / (alpha_d ^ 2));
%             end
%             b = b / max(b);
%             B(:, (ii - 1) * alpha_dnum + jj) = b;
%         end
%     end
%     for ii = 1 : M
%         for m = 1 : M
%             a(m, 1) = exp(1i * pi * (m - 1) * sin(-pi * 0.5 + pi * ii / M));
%         end
%         for jj = 1 : Nt
%             z(ii, 1) = z(ii, 1) + abs(a' * Y(:, jj)) ^ 2;
%         end
%     end
%     z = z - min(z);
%     z = z / max(z);
%     for ii = 1 : alpha_num
%         for jj = 1 : alpha_dnum
%             for m = 1 : M
%                 Lambda((ii - 1) * alpha_dnum + jj) = Lambda((ii - 1) * alpha_dnum + jj) + abs(B(m, (ii - 1) * alpha_dnum + jj) - z(m, 1)) ^ 2;
%             end
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     [~, IX] = sort(Lambda);
%     Eta(1, 1) =  alpha_min + floor((IX(1, 1) - 1) / alpha_dnum)  * preci;
%     temp = rem(IX(1, 1), alpha_dnum);
%     if 0 == temp
%         Eta(2, 1) = alpha_dnum * preci;
%     else
%         Eta(2, 1) = temp * preci;
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     [~, I] = sort(Eta(1, :));
%     eta = Eta(:, I);
else
    Eta = zeros(2, source_num);
    z = zeros(M, 1);
    a = zeros(M, 1);
    b = zeros(M, 1);
    alpha_dnum = theta_d_max / preci;%doa spread range 1`10 degrees
    Lambda = zeros(alpha_dnum, 1);%store of the eigenvalues
    B = zeros(M, alpha_dnum);
    for ii = 1 : M
        for m = 1 : M
            a(m, 1) = exp(1i * pi * (m - 1) * sin(-pi * 0.5 + pi * ii / M));
        end
        for jj = 1 : Nt
            z(ii, 1) = z(ii, 1) + abs(a' * Y(:, jj)) ^ 2;
        end
    end
    z = z - mean(z(1 : 40, 1));
    z = z / max(z);
    [~, I] = max(z);
    Eta(1, 1) = I * 180 / M  - 90;
    for jj = 1 : alpha_dnum
        alpha_d = jj * preci / 180 * pi;
        for m = 1 : M
            theta_m = -pi * 0.5 + pi * m / M;
            b(m, 1) = 1 / sqrt(2 * pi) / alpha_d * exp(-0.5 * (theta_m - Eta(1, 1) / 180 * pi) ^ 2 / (alpha_d ^ 2));
        end
        b = b / max(b);
        B(:, jj) = b;
    end
    for jj = 1 : alpha_dnum
        for m = 1 : M
            Lambda(jj) = Lambda(jj) + abs(B(m, jj) - z(m, 1)) ^ 2;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [~, IX] = min(Lambda);
    Eta(2, 1) = IX * preci;
    %%%%%%%%%%%%%%%%%%%%%%%%%
    c = zeros(M, 1);
    for n = 1 : source_num - 1
        for ii = 1 : M
            theta_m = -pi * 0.5 + pi * ii / M;
            c(ii, 1) = 1 / sqrt(2 * pi) / (Eta(2, n) / 180 * pi) * exp(-0.5 * (theta_m - Eta(1, n) / 180 * pi) ^ 2 / ((Eta(2, n) / 180 * pi) ^ 2));
        end
        z = z - c / max(c);
        z = z / max(z);
        [~, I] = max(z);
        Eta(1, n + 1) = I * 180 / M  - 90;
        for jj = 1 : alpha_dnum
            alpha_d = jj * preci / 180 * pi;
            for m = 1 : M
                theta_m = -pi * 0.5 + pi * m / M;
                b(m, 1) = 1 / sqrt(2 * pi) / alpha_d * exp(-0.5 * (theta_m - Eta(1, n + 1) / 180 * pi) ^ 2 / (alpha_d ^ 2));
            end
            b = b / max(b);
            B(:, jj) = b;
        end
        for jj = 1 : alpha_dnum
            for m = 1 : M
                Lambda(jj) = Lambda(jj) + abs(B(m, jj) - z(m, 1)) ^ 2;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [~, IX] = min(Lambda);
        Eta(2, n + 1) = IX * preci;
    end
    [~, I] = sort(Eta(1, :));
    eta = Eta(:, I);
end


