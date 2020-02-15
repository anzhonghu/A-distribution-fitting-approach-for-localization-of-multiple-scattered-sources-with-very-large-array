function doa_estimate = doasearch_2s(alpha_min, alpha_dnum, preci, IX)
theta_e = zeros(2, 1);
theta_e(1, 1) = alpha_min + floor((IX(1, 1) - 1) / alpha_dnum)  * preci;
temp = rem(IX(1, 1), alpha_dnum);
if 0 == temp
    theta_e(2, 1) = alpha_dnum * preci;
else
    theta_e(2, 1) = temp * preci;
end
doa_estimate = theta_e;