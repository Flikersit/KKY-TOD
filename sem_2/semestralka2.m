clc;
close all;
clear all;
%% inicialization

t = -pi:0.001:pi;
u = [cos(t); sin(t)];
K = 3;

T  = 1;
F  = [1 T; 0 1];           
H  = [1 0];                
Q  = (1/10) * [T^3/3  T^2/2;  
               T^2/2  T   ];
R  = 1;                    
Px0 = [1 1; 1 4];          

x0_real = [5;
           4];

mu_x0_pos = sym('x0_poloha','real');   
mu_x0_vel = sym('x0_rychlost','real');
mu_x0 = [mu_x0_pos; mu_x0_vel];

syms z0 z1

Z = [z0;
    z1];

%% ML
G = [H; 
    H*F];

cov_e_ml = [R 0;
       0 H*Q*H'+R];

inv(cov_e_ml);

x_ml = inv(G' * inv(cov_e_ml)* G)*G'*inv(cov_e_ml)*Z;
x_ml = vpa(x_ml, 4);

estimation_ml = [z0;
                 z1 - z0];

cov_er_ml = inv(G'*inv(cov_e_ml)*G);
   
%% A LMSE
mza = H*x0_real;
P_x0za = Px0*H';
P_za = H*Px0*H' + R;

cov_err_a = Px0 - P_x0za*inv(P_za)*P_x0za';

estimation_a = x0_real + P_x0za*inv(P_za)*(z0 - mza); 

%% B LMSE
mzb = H*F*x0_real;     
P_x0zb = Px0*F'*H';
P_zb = H*F*Px0*F'*H' + H*Q*H' + R;

cov_err_b = Px0 - P_x0zb*inv(P_zb)*P_x0zb';

estimation_b = x0_real + P_x0zb*inv(P_zb)*(z1 - mzb);

%% C LMSE
P_zzc = H*(Px0 + mu_x0*mu_x0')*F'*H' - mu_x0_pos*(mu_x0_pos + mu_x0_vel);
P_zzc = simplify(P_zzc);
P_zc = [P_za P_zzc;
        P_zzc P_zb];
P_zc = double(P_zc);
P_x0zc = [P_x0za P_x0zb];
mzc = [mza; mzb];

inv(vpa(P_zc, 4));
cov_err_c = vpa(Px0 - P_x0zc*inv(P_zc)*P_x0zc', 4);

estimation_c = x0_real + P_x0zc*inv(P_zc)*(Z - mzc);

%% C A 

%bias of the estimator --

R_a = 1/4;

P_z2aa = H*Px0*H' + R_a;
P_z2ba = H*F*Px0*F'*H' + H*Q*H' + R_a;

P_zc2a = [P_z2aa P_zzc;
        P_zzc P_z2ba];
P_zc2a = double(P_zc2a);

inv(P_zc2a)


estimation_ca = simplify(x0_real + P_x0zc*inv(P_zc2a)*(Z - mzc))

%fake error covariance
cov_err_ca_fake = Px0 - P_x0zc*inv(P_zc2a)*P_x0zc';

%real error covariance
cov_err_ca_real = Px0 + P_x0zc*inv(P_zc2a)*P_zc*inv(P_zc2a)*P_x0zc' - P_x0zc*inv(P_zc2a)*P_x0zc' - P_x0zc*inv(P_zc2a)*P_x0zc'


%% C B

% bias of the estimator --

R_b = 4;

P_z2ab = H*Px0*H' + R_b;
P_z2bb = H*F*Px0*F'*H' + H*Q*H' + R_b;

P_zc2b = [P_z2ab P_zzc;
        P_zzc P_z2bb];
P_zc2b = double(P_zc2b);

inv(P_zc2b)

estimation_cb = x0_real + P_x0zc*inv(P_zc2b)*(Z - mzc)

%fake error covariance
cov_err_cb_fake = Px0 - P_x0zc*inv(P_zc2b)*P_x0zc';

%real error covariance
cov_err_cb_real = Px0 + P_x0zc*inv(P_zc2b)*P_zc*inv(P_zc2b)*P_x0zc' - P_x0zc*inv(P_zc2b)*P_x0zc' - P_x0zc*inv(P_zc2b)*P_x0zc'



%% C C

%nestrannost -

mu_x0 = x0_real; 
mu_x0_pos = mu_x0(1);
mu_x0_vel = mu_x0(2);

% change of mz, mx
mzcc = [(mu_x0_pos - 5) + 0*mu_x0_vel;
        (mu_x0_pos + mu_x0_vel -5) + 0*mu_x0_vel];
mxcc = [(mu_x0_pos - 5) + 0*mu_x0_vel;
        mu_x0_vel + 0*mu_x0_vel];

strannost_cc = (mu_x0 - mxcc) - P_x0zc*inv(P_zc)*(mzc - mzcc)
    
estimation_cc = mxcc + P_x0zc*inv(P_zc)*(Z - mzcc)



%% draw elipses a

%fake a
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_ca_fake, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%real a
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_ca_real, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

title('LMSE case A: Measurement noise R = 1/4');
legend('Estimated covariance (model)', 'True covariance (real system)');
xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');


mu_fake = [0; 0];
mu_real = [0; 0];
H2_A = hellinger_distance(mu_fake, cov_err_ca_fake, mu_real, cov_err_ca_real);
fprintf('Case A: Hellinger^2 = %.6f,  H = %.6f\n', H2_A, sqrt(H2_A));

%saveas(gcf, 'elipses_A.eps', 'epsc');
%% draw elipses b

%fake b
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_cb_fake, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%real b
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_cb_real, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

title('LMSE case B: Measurement noise R = 4');
legend('Estimated covariance (model)', 'True covariance (real system)');
xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');


mu_fake = [0; 0];
mu_real = [0; 0];
H2_B = hellinger_distance(mu_fake, cov_err_cb_fake, mu_real, cov_err_cb_real);
fprintf('Case B: Hellinger^2 = %.6f,  H = %.6f\n', H2_B, sqrt(H2_B));

%saveas(gcf, 'elipses_B.eps', 'epsc');
%% draw elipses c

%fake b
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_c, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%real b
x_elipsa = repmat(strannost_cc, 1, length(t)) + K*chol(cov_err_c, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

title('LMSE case C: Mean bias (incorrect state mean)');
legend('Estimated covariance (model)', 'True covariance (real system)');
xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');


mu_fake = [0; 0];
mu_real = strannost_cc;   % Real ellipse is shifted by bias
H2_C = hellinger_distance(mu_fake, cov_err_c, mu_real, cov_err_c);
fprintf('Case C: Hellinger^2 = %.6f,  H = %.6f\n', H2_C, sqrt(H2_C));


%saveas(gcf, 'elipses_C.eps', 'epsc');
%% draw elipses for LMSE + ML


%A
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_a, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%B
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_b, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

%C
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_c, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

%Ac
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_ca_real, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%Bc
x_elipsa = repmat(0 , 1, length(t)) + K*chol(cov_err_cb_real, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

%Cc
x_elipsa = repmat(strannost_cc , 1, length(t)) + K*chol(cov_err_c, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

%Ml
x_elipsa = repmat(0, 1, length(t)) + K*chol(cov_er_ml, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))


title('Comparison of all LMSE and ML estimation cases');

legend('LMSE A (correct)', 'LMSE B (correct)', 'LMSE C (correct)', ...
       'LMSE C (wrong R=1/4)', 'LMSE C (wrong R=4)', ...
       'LMSE C (biased mean)', 'ML estimator');

xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');

%saveas(gcf, 'elipses_all.eps', 'epsc');

%% Compute ellipse areas for all covariance matrices
ellipse_area = @(P, K) pi * K^2 * sqrt(det(P));

A_a  = ellipse_area(cov_err_a,  K);
A_b  = ellipse_area(cov_err_b,  K);
A_c  = ellipse_area(cov_err_c,  K);
A_ca = ellipse_area(cov_err_ca_real, K);
A_cb = ellipse_area(cov_err_cb_real, K);
A_cc = ellipse_area(cov_err_c,  K);     % biased mean uses same covariance
A_ml = ellipse_area(cov_er_ml,  K);

fprintf('\n=== Ellipse areas (scaled by K=%.1f) ===\n', K);
fprintf('LMSE A (correct):      %.4f\n', A_a);
fprintf('LMSE B (correct):      %.4f\n', A_b);
fprintf('LMSE C (correct):      %.4f\n', A_c);
fprintf('LMSE C (wrong R=1/4):  %.4f\n', A_ca);
fprintf('LMSE C (wrong R=4):    %.4f\n', A_cb);
fprintf('LMSE C (biased mean):  %.4f\n', A_cc);
fprintf('ML estimator:          %.4f\n', A_ml);



%% experemential estimation

z_real = zeros(1000, 2);
estimations = zeros(1000, 2, 7);
x0_realization = zeros(1000, 2);

for i = 1:1000
    
    x0_produce = x0_real + chol(Px0, 'lower') * randn(2,1);
    z0_produce = H*x0_produce + normrnd(0, R);
    x1_produce = F*x0_produce + chol(Q, 'lower') * randn(2,1);
    z1_produce = H*x1_produce + normrnd(0, R);
    
    z_real(i, 1) = z0_produce;
    z_real(i, 2) = z1_produce;
    
    x0_realization(i, :) = x0_produce;
    
    x0_A = vpa(subs(estimation_a, {z0}, {z0_produce}), 4);
    x0_B = vpa(subs(estimation_b, {z1}, {z1_produce}), 4);
    x0_C = vpa(subs(estimation_c, {z0, z1}, {z0_produce, z1_produce}), 4);
    
    x0_ML = vpa(subs(estimation_ml, {z0, z1}, {z0_produce, z1_produce}), 4);
    
    x0_Ca = vpa(subs(estimation_ca, {z0, z1}, {z0_produce, z1_produce}), 4);
    x0_Cb = vpa(subs(estimation_cb, {z0, z1}, {z0_produce, z1_produce}), 4);
    x0_Cc = vpa(subs(estimation_cc, {z0, z1}, {z0_produce, z1_produce}), 4);
    
    
    estimations(i, :, 1) = double(x0_A);
    estimations(i, :, 2) = double(x0_B);
    estimations(i, :, 3) = double(x0_C);
    estimations(i, :, 4) = double(x0_Ca);
    estimations(i, :, 5) = double(x0_Cb);
    estimations(i, :, 6) = double(x0_Cc);
    estimations(i, :, 7) = double(x0_ML);
end


%% vel errors

errors_vel = zeros(1000, 7);
for i = 1:7
    errors_vel(:, i) = x0_realization(:, 2) - estimations(:, 2, i);
end


%% histogramm LMSE A

hold on
grid on
h = histogram(errors_vel(:, 1), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 1));
var_emp = var(errors_vel(:, 1), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_a(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE A)');

saveas(gcf, 'histogram_lmsea.eps', 'epsc');

%% histogramm LMSE B

hold on
grid on
h = histogram(errors_vel(:, 2), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 2));
var_emp = var(errors_vel(:, 2), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_b(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE B)');

saveas(gcf, 'histogram_lmseb.eps', 'epsc');

%% histogramm LMSE C

hold on
grid on
h = histogram(errors_vel(:, 3), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 3));
var_emp = var(errors_vel(:, 3), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_c(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE C)');

saveas(gcf, 'histogram_lmsec.eps', 'epsc');

%% histogram Ca


hold on
grid on
h = histogram(errors_vel(:, 4), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 4));
var_emp = var(errors_vel(:, 4), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_ca_fake(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE Cc)');

saveas(gcf, 'histogram_lmseca.eps', 'epsc');

%% histogram Cb

hold on
grid on
h = histogram(errors_vel(:, 5), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 5));
var_emp = var(errors_vel(:, 5), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_cb_fake(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE Cb)');

saveas(gcf, 'histogram_lmsecb.eps', 'epsc');

%% histogram Cc

hold on
grid on
h = histogram(errors_vel(:, 6), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 6));
var_emp = var(errors_vel(:, 6), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_err_c(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (LMSE Cc)');

saveas(gcf, 'histogram_lmsecc.eps', 'epsc');



%% histogram ML

hold on
grid on
h = histogram(errors_vel(:, 7), 40, 'Normalization','pdf');

mu = mean(errors_vel(:, 7));
var_emp = var(errors_vel(:, 7), 0);

x = [-6:.1:8];
y = normpdf(x,mu,sqrt(var_emp));
plot(x, y, 'r', 'LineWidth', 2)


z = normpdf(x, 0, sqrt(cov_er_ml(2, 2)));
plot(x ,z, '--b', 'LineWidth', 2)


legend('Error histogram', ...
       'Empirical Gaussian (method of moments)', ...
       'Theoretical Gaussian');

xlabel('Velocity error');
ylabel('Density');
title('Comparison of empirical and theoretical error densities (ML)');

saveas(gcf, 'histogram_ml.eps', 'epsc');
