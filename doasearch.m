function doa_estimate = doasearch(alpha_min, alpha_dnum, preci, source_num, IX)
theta_e = zeros(2, source_num);
for ii = 1 : source_num
    theta_e(1, ii) = alpha_min + floor((IX(ii, 1) - 1) / alpha_dnum)  * preci;
    temp = rem(IX(ii, 1), alpha_dnum);
    if 0 == temp
        theta_e(2, ii) = alpha_dnum * preci;
    else
        theta_e(2, ii) = temp * preci;
    end
end
[~, I] = sort(theta_e(1, :));
doa_estimate = theta_e(:, I);