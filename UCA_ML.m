clear
array_num = 9;
R = 1;       % Radius is 0.5m
boresight = [45; 45];        % [azimuth; elevation]

array_obj = phased.UCA(array_num, R);
array_azm = (-(array_num - 1)/2 + (1:array_num)' - 1)*360/array_num;        % Elements azimuth
sv_obj = phased.SteeringVector('SensorArray', array_obj);
fc = 300e6;       % f_c = 300MHz
lambda = 3e8/fc;
SteeringVector = @(ang) exp((1j*2*pi*R/lambda)* ...
    (cosd(array_azm)*sind(ang(2))*cosd(ang(1)) + sind(array_azm)*sind(ang(2))*sind(ang(1))));
% sv_0 = sv_obj(fc, boresight);
sv_0 = SteeringVector(boresight);
epsilon_0 = [-sind(boresight(1)), cosd(boresight(2))*cosd(boresight(1));
                cosd(boresight(1)), cosd(boresight(2))*sind(boresight(1));
                0, -sind(boresight(2))];

alpha = array_azm;
alpha(alpha <= 0) = alpha(alpha <= 0) + 90;
alpha(alpha > 0) = alpha(alpha > 0) - 90;
beta_h = [cosd(alpha), sind(alpha), zeros(array_num, 1)];
beta_v = [zeros(array_num, 1), zeros(array_num, 1), ones(array_num, 1)];
A_0 = [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*epsilon_0;
dsv_amz = (1j*2*pi*R/lambda)*(sind(array_azm)*sind(boresight(2))*cosd(boresight(1)) - ...
                              cosd(array_azm)*sind(boresight(2))*sind(boresight(1))).*sv_0;
dsv_elv = (1j*2*pi*R/lambda)*(cosd(array_azm)*cosd(boresight(2))*cosd(boresight(1)) + ...
                              sind(array_azm)*cosd(boresight(2))*sind(boresight(1))).*sv_0;
dEpsilon_amz = [-cosd(boresight(1)), -cosd(boresight(2))*sind(boresight(1));
                -sind(boresight(1)), cosd(boresight(2))*cosd(boresight(1));
                0, 0];
dEpsilon_elv = [0, -sind(boresight(2))*cosd(boresight(1));
                0, -sind(boresight(2))*sind(boresight(1));
                0, -cosd(boresight(2))];
d2sv_aa = -(1j*2*pi*R/lambda)*(sind(array_azm)*sind(boresight(2))*sind(boresight(1)) + ...
                               cosd(array_azm)*sind(boresight(2))*cosd(boresight(1))).*sv_0 ...
          + (1j*2*pi*R/lambda)*(sind(array_azm)*sind(boresight(2))*cosd(boresight(1)) - ...
                              cosd(array_azm)*sind(boresight(2))*sind(boresight(1))).*dsv_amz;
d2sv_ee = -(1j*2*pi*R/lambda)*(cosd(array_azm)*sind(boresight(2))*cosd(boresight(1)) + ...
                               sind(array_azm)*sind(boresight(2))*sind(boresight(1))).*sv_0 ...
          + (1j*2*pi*R/lambda)*(cosd(array_azm)*cosd(boresight(2))*cosd(boresight(1)) + ...
                              sind(array_azm)*cosd(boresight(2))*sind(boresight(1))).*dsv_elv;
d2sv_ae = (1j*2*pi*R/lambda)*(sind(array_azm)*cosd(boresight(2))*cosd(boresight(1)) - ...
                               cosd(array_azm)*cosd(boresight(2))*sind(boresight(1))).*sv_0 ...
          + (1j*2*pi*R/lambda)*(sind(array_azm)*sind(boresight(2))*cosd(boresight(1)) - ...
                              cosd(array_azm)*sind(boresight(2))*sind(boresight(1))).*dsv_elv;
d2Epsilon_aa = [sind(boresight(1)), -cosd(boresight(2))*cosd(boresight(1));
                -cosd(boresight(1)), -cosd(boresight(2))*sind(boresight(1));
                0, 0];
d2Epsilon_ee = [0, -cosd(boresight(2))*cosd(boresight(1));
                0, -cosd(boresight(2))*sind(boresight(1));
                0, sind(boresight(2))];
d2Epsilon_ae = [0, sind(boresight(2))*sind(boresight(1));
                0, -sind(boresight(2))*cosd(boresight(1));
                0, 0];
            
dA_a = [diag(dsv_amz)*beta_h; diag(dsv_amz)*beta_v]*epsilon_0 + ...
       [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*dEpsilon_amz;
dA_e = [diag(dsv_elv)*beta_h; diag(dsv_elv)*beta_v]*epsilon_0 + ...
       [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*dEpsilon_elv;
d2A_aa = [diag(d2sv_aa)*beta_h; diag(d2sv_aa)*beta_v]*epsilon_0 + ...
         [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*d2Epsilon_aa + ...
         2*[diag(dsv_amz)*beta_h; diag(dsv_amz)*beta_v]*dEpsilon_amz;
d2A_ee = [diag(d2sv_ee)*beta_h; diag(d2sv_ee)*beta_v]*epsilon_0 + ...
         [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*d2Epsilon_ee + ...
         2*[diag(dsv_elv)*beta_h; diag(dsv_elv)*beta_v]*dEpsilon_elv;
d2A_ae = [diag(d2sv_ae)*beta_h; diag(d2sv_ae)*beta_v]*epsilon_0 + ...
         [diag(sv_0)*beta_h; diag(sv_0)*beta_v]*d2Epsilon_ae + ...
         [diag(dsv_amz)*beta_h; diag(dsv_amz)*beta_v]*dEpsilon_elv + ...
         [diag(dsv_elv)*beta_h; diag(dsv_elv)*beta_v]*dEpsilon_amz;


JNR = 55;
SNAPSHOTS = 200;
DOA_j = [25; -25];
epsilon_j = [-sind(DOA_j(1)), cosd(DOA_j(2))*cosd(DOA_j(1));
              cosd(DOA_j(1)), cosd(DOA_j(2))*sind(DOA_j(1));
              0, -sind(DOA_j(2))];
% sv_j = sv_obj(fc, DOA_j);
sv_j = SteeringVector(DOA_j);
polar_j = deg2rad([40; 30]);
h_j = [cos(polar_j(1)); sin(polar_j(1))*exp(1j*polar_j(2))];
svv_j = diag(sv_j)*beta_v*epsilon_j*h_j;
svh_j = diag(sv_j)*beta_h*epsilon_j*h_j;
jammer = sqrt(10^(JNR/10)).*exp(1j*2*pi*rand(1, SNAPSHOTS));
noise = sqrt(.5)*randn(array_num, SNAPSHOTS) + 1j*sqrt(.5)*randn(array_num, SNAPSHOTS);
xh_jn = svh_j*jammer + noise;
noise = sqrt(.5)*randn(array_num, SNAPSHOTS) + 1j*sqrt(.5)*randn(array_num, SNAPSHOTS);
xv_jn = svv_j*jammer + noise;
x_jn = [xh_jn; xv_jn];
cov_mat = x_jn*x_jn'/SNAPSHOTS;

W = cov_mat\A_0/sqrtm(A_0'/cov_mat*A_0);
D_a = cov_mat\dA_a/sqrtm(A_0'/cov_mat*A_0);
mu_a = (A_0'/cov_mat*A_0)\real(dA_a'/cov_mat*A_0);
W_a = D_a - W*mu_a;
dD_aa = cov_mat\d2A_aa/sqrtm(A_0'/cov_mat*A_0) - D_a*mu_a;
dmu_aa = -2*pinv(A_0'/cov_mat*A_0)*pinv(A_0'/cov_mat*A_0) * ...
        real(dA_a'/cov_mat*A_0)*real(dA_a'/cov_mat*A_0) + (A_0'/cov_mat*A_0)\ ...
        real(d2A_aa'/cov_mat*A_0 + dA_a'/cov_mat*dA_a);
W_aa = dD_aa - W_a*mu_a - W*dmu_aa;
D_e = cov_mat\dA_e/sqrtm(A_0'/cov_mat*A_0);
mu_e = (A_0'/cov_mat*A_0)\real(dA_e'/cov_mat*A_0);
W_e = D_e - W*mu_e;
dD_ee = cov_mat\d2A_ee/sqrtm(A_0'/cov_mat*A_0) - D_e*mu_e;
dmu_ee = -2*pinv(A_0'/cov_mat*A_0)*pinv(A_0'/cov_mat*A_0) * ...
        real(dA_e'/cov_mat*A_0)*real(dA_e'/cov_mat*A_0) + (A_0'/cov_mat*A_0)\ ...
        real(d2A_ee'/cov_mat*A_0 + dA_e'/cov_mat*dA_e);
W_ee = dD_ee - W_e*mu_e - W*dmu_ee;
dD_ae = cov_mat\d2A_ae/sqrtm(A_0'/cov_mat*A_0) - D_a*mu_e;
dmu_ae = -2*pinv(A_0'/cov_mat*A_0)*pinv(A_0'/cov_mat*A_0) * ...
        real(dA_e'/cov_mat*A_0)*real(dA_a'/cov_mat*A_0) + (A_0'/cov_mat*A_0)\ ...
        real(d2A_ae'/cov_mat*A_0 + dA_a'/cov_mat*dA_e);
W_ae = dD_ae - W_e*mu_a - W*dmu_ae;


DOA_s = boresight + [-3; 3];
% sv_s = sv_obj(fc, DOA_s);
sv_s = SteeringVector(DOA_s);
epsilon_s = [-sind(DOA_s(1)), cosd(DOA_s(2))*cosd(DOA_s(1));
                cosd(DOA_s(1)), cosd(DOA_s(2))*sind(DOA_s(1));
                0, -sind(DOA_s(2))];
gamma = deg2rad(0:90);
eta = 0;
SNR = 15;
MC_L = 1000;
RMSE = zeros(2, length(gamma));
RMSE_p = zeros(2, length(gamma));
for l = 1:MC_L
    for n = 1:length(gamma)
        h_s = [cos(gamma(n)); sin(gamma(n))*exp(1j*eta)];
        svh_s = diag(sv_s)*beta_h*epsilon_s*h_s;
        svv_s = diag(sv_s)*beta_v*epsilon_s*h_s;
        signal = sqrt(10^(SNR/10)).*exp(1j*2*pi*rand(1, SNAPSHOTS));
        jammer = sqrt(10^(JNR/10)).*exp(1j*2*pi*rand(1, SNAPSHOTS));
        noise = sqrt(.5)*randn(array_num, SNAPSHOTS) + 1j*sqrt(.5)*randn(array_num, SNAPSHOTS);
        xh = svh_s*signal + svh_j*jammer + noise;
        noise = sqrt(.5)*randn(array_num, SNAPSHOTS) + 1j*sqrt(.5)*randn(array_num, SNAPSHOTS);
        xv = svv_s*signal + svv_j*jammer + noise;  
        x = [xh; xv];
        
        dir = 0;
        for m = 1:SNAPSHOTS
            dlf_a = 2*real((x(:, m)'*W_a*W'*x(:, m))./(x(:, m)'*W*W'*x(:, m)));
            dlf_e = 2*real((x(:, m)'*W_e*W'*x(:, m))./(x(:, m)'*W*W'*x(:, m)));
            d2lf_aa = 2*real((x(:, m)'*W_aa*W'*x(:, m) + x(:, m)'*W_a*W_a'*x(:, m))./ ...
                             (x(:, m)'*W*W'*x(:, m))) - dlf_a^2;
            d2lf_ee = 2*real((x(:, m)'*W_ee*W'*x(:, m) + x(:, m)'*W_e*W_e'*x(:, m))./ ...
                             (x(:, m)'*W*W'*x(:, m))) - dlf_e^2;
            d2lf_ae = 2*real((x(:, m)'*W_ae*W'*x(:, m) + x(:, m)'*W_a*W_e'*x(:, m))./ ...
                             (x(:, m)'*W*W'*x(:, m))) - dlf_a*dlf_e;
            H = [d2lf_aa, d2lf_ae; d2lf_ae, d2lf_ee];
            dir = dir + rad2deg(H\[dlf_a; dlf_e]);
        end
        DOA_hat = boresight - dir./SNAPSHOTS;
       
        RMSE(:, n) = RMSE(:, n) + (DOA_hat - DOA_s).^2;
        
        sv_hat = SteeringVector(DOA_hat);
        epsilon_hat = [-sind(DOA_hat(1)), cosd(DOA_hat(2))*cosd(DOA_hat(1));
                        cosd(DOA_hat(1)), cosd(DOA_hat(2))*sind(DOA_hat(1));
                        0, -sind(DOA_hat(2))];
        A_hat = [diag(sv_hat)*beta_h; diag(sv_hat)*beta_v]*epsilon_hat;
        e_hat = (A_hat'/cov_mat*A_hat)\A_hat'/cov_mat*x;
        C_e = e_hat*e_hat'/SNAPSHOTS;
        [eigVecs, eigVals] = eig(C_e);
        eigVals = diag(eigVals);
        [~, idx] = max(eigVals);
        h_hat = eigVecs(:, idx);
        eta_hat = angle(h_hat(2)/h_hat(1));
        gamma_hat = atan(real(h_hat(2)*exp(-1j*eta_hat)/h_hat(1)));
        RMSE_p(:, n) = RMSE_p(:, n) + rad2deg([gamma(n); eta] - [gamma_hat; eta_hat]).^2;
    end
end
RMSE = sqrt(RMSE./MC_L);
RMSE_p = sqrt(RMSE_p./MC_L);

figure
plot(rad2deg(gamma), RMSE(1, :))
grid on
xlabel("SNR (dB)")
ylabel("RMSE (degree)")
% ylim([0, 2])
% zlabel("RMSE (degree)")
title("Azimuth RMSE")

figure
plot(rad2deg(gamma), RMSE(2, :))
grid on
xlabel("SNR (dB)")
ylabel("RMSE (degree)")
% ylim([0, 2])
% zlabel("RMSE (degree)")
title("Elevation RMSE")

figure
plot(rad2deg(gamma), RMSE_p(1, :))
grid on
xlabel("\gamma (degree)")
ylabel("RMSE (degree)")
% ylim([0, 2])
% zlabel("RMSE (degree)")
title("\gamma RMSE")

figure
plot(rad2deg(gamma), RMSE_p(2, :))
grid on
xlabel("\gamma (degree)")
ylabel("RMSE (degree)")
% ylim([0, 2])
% zlabel("RMSE (degree)")
title("\eta RMSE")

