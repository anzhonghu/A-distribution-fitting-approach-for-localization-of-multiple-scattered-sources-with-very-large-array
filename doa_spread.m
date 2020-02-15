%This file is for DOA estimation in LS-MIMO angular spread scenarios

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10;%antenna number
Nt = 100;%number of snapshots
preci = 0.1;%resolution of 0.1 degree, preci>=180/M
Y = zeros(M, Nt);%received signal
Sk = zeros(1, Nt);%transmitted signal for each terminal
N = zeros(M, Nt);%received noise at the BS
a = zeros(M, 1);%steering vector
rmse1_bound = zeros(2, 1);
rmse2_bound = zeros(2, 1);
Ite_num = 100;%number of iterations
phi = zeros(2, 4);
theta_d_max = 12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spread_n = 5;
spread = [2;4;6;8;10;];
rmse_store = zeros(2 * spread_n, 7);
SNR = 20;%dB
K = 1;%user number
Nk = 50 * ones(K, 1);%number of paths
theta = 0;%nominal angles, in ascending order
theta_d = 3;%angular deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1 : spread_n
    theta_d = spread(nn, 1);
    rmse = zeros(2, 6);%2 parameters, 6 approaches
    rmse1 = zeros(2, 4);%2 parameters, 4 approaches
    rmse2 = zeros(2, 4);%2 parameters, 4 approaches
    for ii = 1 : Ite_num
        %generate the received signal
        Y = zeros(M, Nt);%received signal
        N = (randn(M, Nt) + 1i * randn(M, Nt)) / sqrt(2);
        Y = Y + N;
        k = 1;
        amp_k = sqrt(10 ^ (SNR * 0.1) / Nk(k, 1));
        Sk = sign(randn(1, Nt));%BPSK modulation
        for jj = 1 : Nt
            if 0 == Sk(1, jj)
                Sk(1, jj) = 1;
            else
            end
            Sk(1, jj) = Sk(1, jj) * amp_k;
            alpha_k = (randn(Nk(k, 1), 1) + 1i * randn(Nk(k, 1), 1)) / sqrt(2);%small scale fading
            theta_k = theta(k, 1) * ones(Nk(k, 1), 1);
            theta_k = theta_k + randn(Nk(k, 1), 1) * theta_d(k, 1);%gaussian distribution
            for kk = 1 : Nk(k, 1)
                for m = 1 : M
                    a(m, 1) = exp(1i * pi * (m - 1) * sin(theta_k(kk, 1) / 180 * pi));
                end
                Y(:, jj) = Y(:, jj) + alpha_k(kk, 1) * a * Sk(1, jj);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %estimation
        RY = zeros(M, M);
        for jj = 1 : Nt
            RY = RY + Y(:, jj) * Y(:, jj)';
        end
        RY = RY / Nt;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DISPARE
        xi = eigsearch(M, theta, RY, preci, theta_d_max);
        rmse(:, 1) = rmse(:, 1) + (xi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
        %GC approach, Subspace approach
        phi = matrixsearch(M, theta, RY, preci, theta_d_max);
        rmse(:, 2) = rmse(:, 2) + (phi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
        rmse(:, 3) = rmse(:, 3) + (phi(:, 2) - [theta(1, 1); theta_d(1, 1)]).^2;
        %Subdiag approach, need only estimate the nominal angle, because
        %the estimation of angular spread is multidimensional
        psi = subdiagsearch(M, theta, RY, preci);
        rmse(1, 4) = rmse(1, 4) + (psi(:, 1) - theta(1, 1))^2;
        %COMET-EXIP approach, need only estimate the nominal angle,  because
        %the estimation of angular spread is multidimensional
        gamma = comet(M, theta, RY, preci);
        rmse(1, 5) = rmse(1, 5) + (gamma(:, 1) -  theta(1, 1))^2;
        %Proposed approach
        eta = vectorsearch(Nt, M, theta, Y, preci, theta_d_max);
        rmse(:, 6) = rmse(:, 6) + (eta(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
    end
    %%%%%%%%%%%%%%%%%%%
    %CRB
    bound = crb(Nt, K, M, theta, theta_d, SNR);
    rmse_store(nn, :) = [sqrt( rmse(1, :) / Ite_num), bound(1, 1);];
    rmse_store(nn + spread_n, :) = [sqrt( rmse(2, :) / Ite_num), bound(2, 1);];
    sprintf('%d', nn)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
plot(spread,rmse_store(1:spread_n,1),'g-.*','LineWidth',2,'MarkerSize',10)
hold on
plot(spread,rmse_store(1:spread_n,2),'m-.o','LineWidth',2,'MarkerSize',14)
plot(spread,rmse_store(1:spread_n,3),'b-.x','LineWidth',2,'MarkerSize',10)
plot(spread,rmse_store(1:spread_n,4),'k-d','LineWidth',2,'MarkerSize',14)
plot(spread,rmse_store(1:spread_n,5),'r-','LineWidth',2,'MarkerSize',14)
plot(spread,rmse_store(1:spread_n,6),'ko-','LineWidth',2,'MarkerSize',14)
plot(spread,rmse_store(1:spread_n,7),'g-','LineWidth',2,'MarkerSize',14)
grid on
%le = legend('Pilot-based estimator','Estimator of [11]','Estimator of [19]','Proposed estimator', 'Location','Southwest');
% set(le,'Fontsize',14,'Fontname','Times')
xlabel('DOA spread','Fontsize',16,'Fontname','Times')
ylabel('RMSE','Fontsize',20,'Fontname','Times')
%print(h,'-dpdf','BER_SNR_K')




