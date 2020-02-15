%This file is for DOA estimation in LS-MIMO angular spread scenarios

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 100;%number of snapshots
preci = 0.1;%resolution of 0.1 degree, preci>=180/M
theta_range = 4;%the range of searched angle, degree
Sk = zeros(1, Nt);%transmitted signal for each terminal
rmse1_bound = zeros(2, 1);
rmse2_bound = zeros(2, 1);
Ite_num = 1e2;%number of iterations
phi = zeros(2, 4);
theta_d_max = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 2;%user number
SNR = 10;%dB
Nk = 50 * ones(K, 1);%number of paths
if 1 == K
    M_n = 5;
    Ms = [2;4;6;8;10];%dB
    rmse_store = zeros(2 * M_n, 7);
    theta = 0;%nominal angles, in ascending order
    theta_d = 3;%angular deviation
else
    M_n = 8;
    Ms = [6:2:20]';%dB
    rmse_store = zeros(2 * M_n, 1);
    theta = [-4; 9];%nominal angles, in ascending order
    theta_d = [2; 3];%angular deviation
    bound = zeros(4, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1 : M_n
    M = Ms(nn, 1);
    Y = zeros(M, Nt);%received signal
    a = zeros(M, 1);%steering vector
    N = zeros(M, Nt);%received noise at the BS
    bound = crb(Nt, K, M, theta, theta_d, SNR);
    rmse_store(nn, :) = 0.5 * (bound(1, 1) + bound(2, 1));
    rmse_store(nn + M_n, :) = 0.5 * (bound(3, 1) + bound(4, 1));
end




