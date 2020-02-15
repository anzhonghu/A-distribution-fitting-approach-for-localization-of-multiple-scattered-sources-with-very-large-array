%This file is for DOA estimation in LS-MIMO angular spread scenarios

clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt = 300;%number of snapshots
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
SNR = 20;%dB
Nk = 50 * ones(K, 1);%number of paths
if 1 == K
    M_n = 5;
    Ms = [50;100;150;200;250];%dB
    rmse_store = zeros(2 * M_n, 7);
    theta = 0;%nominal angles, in ascending order
    theta_d = 3;%angular deviation
else
    M_n = 5;
    Ms = [50;100;150;200;250];%dB
    rmse_store = zeros(2 * M_n, 5);
    rmse_store1 = zeros(2 * M_n, 5);
    rmse_store2 = zeros(2 * M_n, 5);
    theta = [0; 20];%nominal angles, in ascending order
    theta_d = [2; 3];%angular deviation
    bound = zeros(4, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn = 1 : M_n
    M = Ms(nn, 1);
    Y = zeros(M, Nt);%received signal
    a = zeros(M, 1);%steering vector
    N = zeros(M, Nt);%received noise at the BS
    rmse = zeros(2, 6);%2 parameters, 6 approaches
    rmse1 = zeros(2, 4);%2 parameters, 4 approaches
    rmse2 = zeros(2, 4);%2 parameters, 4 approaches
%     for ii = 1 : Ite_num
%         %generate the received signal
%         Y = zeros(M, Nt);%received signal
%         N = (randn(M, Nt) + 1i * randn(M, Nt)) / sqrt(2);
%         Y = Y + N;
%         for k = 1 : K
%             amp_k = sqrt(10 ^ (SNR * 0.1) / Nk(k, 1));
%             Sk = sign(randn(1, Nt));%BPSK modulation
%             for jj = 1 : Nt
%                 if 0 == Sk(1, jj)
%                     Sk(1, jj) = 1;
%                 else
%                 end
%                 Sk(1, jj) = Sk(1, jj) * amp_k;
%                 alpha_k = (randn(Nk(k, 1), 1) + 1i * randn(Nk(k, 1), 1)) / sqrt(2);%small scale fading
%                 theta_k = theta(k, 1) * ones(Nk(k, 1), 1);
%                 theta_k = theta_k + randn(Nk(k, 1), 1) * theta_d(k, 1);%gaussian distribution
%                 for kk = 1 : Nk(k, 1)
%                     for m = 1 : M
%                         a(m, 1) = exp(1i * pi * (m - 1) * sin(theta_k(kk, 1) / 180 * pi));
%                     end
%                     Y(:, jj) = Y(:, jj) + alpha_k(kk, 1) * a * Sk(1, jj);
%                 end
%             end
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %estimation
%         RY = zeros(M, M);
%         for jj = 1 : Nt
%             RY = RY + Y(:, jj) * Y(:, jj)';
%         end
%         RY = RY / Nt;
%         %%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if 1 == K
%             %DISPARE
% %             xi = eigsearch(M, theta, RY, preci, theta_d_max, theta_range);
% %             rmse(:, 1) = rmse(:, 1) + (xi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
% %             %GC approach, Subspace approach
% %             phi = matrixsearch(M, theta, RY, preci, theta_d_max, theta_range);
% %             rmse(:, 2) = rmse(:, 2) + (phi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
% %             rmse(:, 3) = rmse(:, 3) + (phi(:, 2) - [theta(1, 1); theta_d(1, 1)]).^2;
% %             %Subdiag approach, need only estimate the nominal angle, because
% %             %the estimation of angular spread is multidimensional
% %             psi = subdiagsearch(M, theta, RY, preci, theta_range);
% %             rmse(1, 4) = rmse(1, 4) + (psi(:, 1) - theta(1, 1))^2;
% %             %COMET-EXIP approach, need only estimate the nominal angle,  because
% %             %the estimation of angular spread is multidimensional
% %             gamma = comet(M, theta, RY, preci, theta_range);
% %             rmse(1, 5) = rmse(1, 5) + (gamma(:, 1) -  theta(1, 1))^2;
%             %Proposed approach
%             eta = vectorsearch_n(Nt, M, theta, Y, preci, theta_d_max, theta_range);
%             rmse(:, 6) = rmse(:, 6) + (eta(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
%         else
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %DISPARE
%             xi = eigsearch(M, theta, RY, preci, theta_d_max, theta_range);
%             rmse1(:, 1) = rmse1(:, 1) + (xi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
%             %GC approach, Subspace approach
%             phi = matrixsearch(M, theta, RY, preci, theta_d_max, theta_range);
%             rmse1(:, 2) = rmse1(:, 2) + (phi(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
%             rmse1(:, 3) = rmse1(:, 3) + (phi(:, 3) - [theta(1, 1); theta_d(1, 1)]).^2;
%             %Proposed approach
%             eta = vectorsearch_n(Nt, M, theta, Y, preci, theta_d_max, theta_range);
%             rmse1(:, 4) = rmse1(:, 4) + (eta(:, 1) - [theta(1, 1); theta_d(1, 1)]).^2;
%             %%%%%%%%%%%%%%%%%%%%%%%%%
%             %DISPARE
%             rmse2(:, 1) = rmse2(:, 1) + (xi(:, 2) - [theta(2, 1); theta_d(2, 1)]).^2;
%             %GC approach, Subspace approach
%             rmse2(:, 2) = rmse2(:, 2) + (phi(:, 2) - [theta(2, 1); theta_d(2, 1)]).^2;
%             rmse2(:, 3) = rmse2(:, 3) + (phi(:, 4) - [theta(2, 1); theta_d(2, 1)]).^2;
%             %Proposed approach
%             rmse2(:, 4) = rmse2(:, 4) + (eta(:, 2) - [theta(2, 1); theta_d(2, 1)]).^2;
%         end
%     end
    %%%%%%%%%%%%%%%%%%%
    %CRB
    bound = crb(Nt, K, M, theta, theta_d, SNR);
    if 1 == K
        rmse_store(nn, :) = [sqrt(rmse(1, :) / Ite_num), bound(1, 1);];
        rmse_store(nn + M_n, :) = [sqrt(rmse(2, :) / Ite_num), bound(2, 1);];
    else
        rmse_store1(nn, :) = [sqrt(rmse1(1, :) / Ite_num), bound(1, 1)];
        rmse_store1(nn + M_n, :) = [sqrt(rmse1(2, :) / Ite_num), bound(3, 1);];
        rmse_store2(nn, :) = [sqrt(rmse2(1, :) / Ite_num), bound(2, 1)];
        rmse_store2(nn + M_n, :) = [sqrt(rmse2(2, :) / Ite_num), bound(4, 1);];
        rmse_store(nn, :) = (rmse_store1(nn, :) + rmse_store2(nn, :)) * 0.5;
        rmse_store(nn + M_n, :) = (rmse_store1(nn + M_n, :) + rmse_store2(nn + M_n, :)) * 0.5;
    end
    sprintf('%d', nn)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
xx = axes('FontSize',16);
if 1 == K
    plot(Ms,rmse_store(1:M_n,1),'g-.*','LineWidth',2,'MarkerSize',10)
    hold on
    plot(Ms,rmse_store(1:M_n,2),'m-.o','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(1:M_n,3),'b-.x','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(1:M_n,4),'k-d','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(1:M_n,5),'r-','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(1:M_n,6),'ko-','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(1:M_n,7),'g-','LineWidth',2,'MarkerSize',14)
else
    plot(Ms,rmse_store(1:M_n,1),'g-.o','LineWidth',2,'MarkerSize',14)
    hold on
    plot(Ms,rmse_store(1:M_n,2),'m-.x','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(1:M_n,3),'b-.*','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(1:M_n,4),'k-d','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(1:M_n,5),'r-','LineWidth',2,'MarkerSize',14)
end
grid on
le = legend('DISPARE','GC estimator','Subspace estimator','Proposed estimator', 'CRB', 'Location','Northwest');
% le = legend('散射信号参数估计方法','波束成形方法','子空间方法','本发明方法', '克拉美罗界', 'Location','Northwest');
set(le,'Fontsize',14,'Fontname','Times')
ylim([0, 0.4])
xlabel('Number of BS antennas','Fontsize',16,'Fontname','Times')
ylabel('RMSE (degrees)','Fontsize',20,'Fontname','Times')
% xlabel('基站天线数量','Fontsize',16,'Fontname','Times')
% ylabel('RMSE (度)','Fontsize',20,'Fontname','Times')
%print(h,'-dpdf','RMSE_DOA_antenna')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h1 = figure;
set(h1,'PaperType','A4');
yy = axes('FontSize',16);
if 1 == K
    plot(Ms,rmse_store(M_n+1:2*M_n,1),'g-.*','LineWidth',2,'MarkerSize',10)
    hold on
    plot(Ms,rmse_store(M_n+1:2*M_n,2),'m-.o','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(M_n+1:2*M_n,3),'b-.x','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(M_n+1:2*M_n,4),'k-d','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(M_n+1:2*M_n,5),'r-','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(M_n+1:2*M_n,6),'ko-','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(M_n+1:2*M_n,7),'g-','LineWidth',2,'MarkerSize',14)
else
    plot(Ms,rmse_store(M_n+1:2*M_n,1),'g-.o','LineWidth',2,'MarkerSize',14)
    hold on
    plot(Ms,rmse_store(M_n+1:2*M_n,2),'m-.x','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(M_n+1:2*M_n,3),'b-.*','LineWidth',2,'MarkerSize',10)
    plot(Ms,rmse_store(M_n+1:2*M_n,4),'k-d','LineWidth',2,'MarkerSize',14)
    plot(Ms,rmse_store(M_n+1:2*M_n,5),'r-','LineWidth',2,'MarkerSize',14)
end
grid on
le = legend('DISPARE','GC estimator','Subspace estimator','Proposed estimator', 'CRB', 'Location','Northwest');
% le = legend('散射信号参数估计方法','波束成形方法','子空间方法','本发明方法', '克拉美罗界', 'Location','Northwest');
set(le,'Fontsize',14,'Fontname','Times')
ylim([0, 1])
xlabel('Number of BS antennas','Fontsize',16,'Fontname','Times')
ylabel('RMSE (degrees)','Fontsize',20,'Fontname','Times')
% xlabel('基站天线数量','Fontsize',16,'Fontname','Times')
% ylabel('RMSE (度)','Fontsize',20,'Fontname','Times')
%print(h1,'-dpdf','RMSE_spread_antenna')



