clc;
clear all;
close all;

%% inicialization

syms z0 z1 z2 z3 

t = -pi:0.001:pi;
u = [cos(t); sin(t)];

m = [5; 4];
K = 3;


T = 1;
R = 1;
q = 0.1;

H = [1 0];
F = [1 T;
     0 1];
%dotaz
J = [H*F^0;
     H*F^1;
     H*F^2;
     H*F^3];

cov_e = [R 0 0 0;
    0 q*((T^3)/3)+R q*((5*(T^3))/6) q*((4*(T^3))/3);
    0 q*((5*(T^3))/6) q*(8*(T^3)/3)+R q*((14*(T^3))/3);
    0 q*((4*(T^3))/3) q*((14*(T^3))/3) q*14*(T^3) + R];

%% point F
W_0123 = inv(cov_e)

J = [H*F^0;
     H*F^1;
     H*F^2;
     H*F^3];

Z = [z0;
     z1;
     z2;
     z3];

x0_0123 = inv(J'*W_0123*J)* J' * W_0123 * Z
Cov_er_F = inv(J'*W_0123*J);

%% point A
Cov_e_01 = [cov_e(1, 1) cov_e(1, 2);
            cov_e(2, 1) cov_e(2, 2)];
W_01 = inv(Cov_e_01)

J = [H*F^0;
     H*F^1];
 
Z = [z0;
     z1];
 
x0_01 = inv(J'*W_01*J)* J' * W_01 * Z;
vpa(x0_01, 4)

Cov_er_A = inv(J'*W_01*J);

%% point B 

Cov_e_02 = [cov_e(1, 1) cov_e(1, 3);
            cov_e(3, 1) cov_e(3, 3)];
W_02 = inv(Cov_e_02)

J = [H*F^0;
     H*F^2];
 
Z = [z0;
     z2];
 
x0_02 = inv(J'*W_02*J)* J' * W_02 * Z;
vpa(x0_02, 4)
Cov_er_B = inv(J'*W_02*J);
%% point C

Cov_e_12 = [cov_e(2, 2) cov_e(2, 3);
            cov_e(3, 2) cov_e(3, 3)];
W_12 = inv(Cov_e_12)

J = [H*F^1;
     H*F^2];
 
Z = [z1;
     z2];
 
x0_12 = inv(J'*W_12*J)* J' * W_12 * Z;
vpa(x0_12, 4)
Cov_er_C = inv(J'*W_12*J);
%% point D
Cov_e_012 = [cov_e(1, 1) cov_e(2, 1) cov_e(3, 1);
             cov_e(1, 2) cov_e(2, 2) cov_e(3, 2);
             cov_e(1, 3) cov_e(2, 3) cov_e(3, 3)];
W_012 = inv(Cov_e_012)

J = [H*F^0;
     H*F^1;
     H*F^2];
 
Z = [z0;
     z1;
     z2];
 
x0_012 = inv(J'*W_012*J)* J' * W_012 * Z;
vpa(x0_012, 4)
Cov_er_D = inv(J'*W_012*J);
%% point E
Cov_e_123 = [cov_e(2, 2) cov_e(2, 3) cov_e(2, 4);
             cov_e(3, 2) cov_e(3, 3) cov_e(3, 4);
             cov_e(4, 2) cov_e(4, 3) cov_e(4, 4)];
W_123 = inv(Cov_e_123)

J = [H*F^1;
     H*F^2;
     H*F^3];
 
Z = [z1;
     z2;
     z3];
 
x0_123 = inv(J'*W_123*J)* J' * W_123 * Z;
vpa(x0_123, 4)
Cov_er_E = inv(J'*W_123*J);
%% ellipse 
%A

x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_A, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))


%B
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_B, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))
%C
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_C, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))
%D
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_D, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))
%E 
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_E, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))
%F
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_F, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))


xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');
title('Theoretical error covariance ellipses');
legend('A', 'B', 'C', 'D', 'E', 'F');

saveas(gcf, 'elipses.eps', 'epsc');

%% realization of stochastic process 

Q = q * [1/3 1/2; 1/2 1];
z_real = zeros(1000, 4);
estimation = zeros(1000, 2, 6);

for i = 1:1000
    z_0 = H*m + normrnd(0, R);
    x_1 = F*m + chol(Q, 'lower') * randn(2,1);
    z_1 = H*x_1 + normrnd(0, R);
    x_2 = F*x_1 + chol(Q, 'lower') * randn(2,1);
    z_2 = H*x_2 + normrnd(0, R);
    x_3 = F*x_2 + chol(Q, 'lower') * randn(2,1);
    z_3 = H*x_3 + normrnd(0, R);
    
    z_real(i, 1) = z_0;
    z_real(i, 2) = z_1;
    z_real(i, 3) = z_2;
    z_real(i, 4) = z_3;
    
    x0_A = vpa(subs(x0_01, {z0, z1}, {z_0, z_1}), 4);
    x0_B = vpa(subs(x0_02, {z0, z2}, {z_0, z_2}), 4);
    x0_C = vpa(subs(x0_12, {z1, z2}, {z_1, z_2}), 4);
    x0_D = vpa(subs(x0_012, {z0, z1, z2}, {z_0, z_1, z_2}), 4);
    x0_E = vpa(subs(x0_123, {z1, z2, z3}, {z_1, z_2, z_3}), 4);
    x0_F = vpa(subs(x0_0123, {z0, z1, z2, z3}, {z_0, z_1, z_2, z_3}), 4);
    
    estimation(i, :, 1) = double(x0_A);
    estimation(i, :, 2) = double(x0_B);
    estimation(i, :, 3) = double(x0_C);
    estimation(i, :, 4) = double(x0_D);
    estimation(i, :, 5) = double(x0_E);
    estimation(i, :, 6) = double(x0_F);

end


%% histogram D
grid on
hold on
h = histogram(estimation(:, 2, 4), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 4));
sigma = sqrt(Cov_er_D(2, 2));


x = [1:.1:7];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point D');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_D.eps', 'epsc');
%% histogram A
grid on
hold on
h = histogram(estimation(:, 2, 1), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 1));
sigma = sqrt(Cov_er_A(2, 2));


x = [0:.1:7.5];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point A');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_A.eps', 'epsc');
%% histogram B
grid on
hold on
h = histogram(estimation(:, 2, 2), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 2));
sigma = sqrt(Cov_er_B(2, 2));


x = [0:.1:7.5];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point B');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_B.eps', 'epsc');
%% histogram C
grid on
hold on
h = histogram(estimation(:, 2, 3), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 3));
sigma = sqrt(Cov_er_C(2, 2));


x = [0:.1:7.5];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point C');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_C.eps', 'epsc');
%% histogram E
grid on
hold on
h = histogram(estimation(:, 2, 5), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 5));
sigma = sqrt(Cov_er_E(2, 2));


x = [0:.1:7.5];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point E');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_E.eps', 'epsc');
%% histogram F
grid on
hold on
h = histogram(estimation(:, 2, 6), 40, 'Normalization','pdf')

mu = mean(estimation(:, 2, 6));
sigma = sqrt(Cov_er_F(2, 2));

x = [0:.1:7.5];
y = normpdf(x,mu,sigma);
plot(x, y)

xlabel('Estimated x_2 (velocity)');
ylabel('Probability density');
title('Histogram of estimated x_2 for point F');
legend('Empirical distribution', 'Theoretical normal PDF');

saveas(gcf, 'histogram_F.eps', 'epsc');
%% ellipse 
X_hat_F = squeeze(estimation(:, :, 6));
E_F = X_hat_F - m';              
Cov_est_F = cov(E_F);

%estimation
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_est_F, 'lower') * u;
grid on;
hold on;
plot(x_elipsa(1, :), x_elipsa(2, :))

%theoretical
x_elipsa = repmat(0 , 1, length(t)) + K*chol(Cov_er_F, 'lower') * u;
plot(x_elipsa(1, :), x_elipsa(2, :))

legend('estimation', 'theoretical');
xlabel('e_1 (position error)');
ylabel('e_2 (velocity error)');
title('Estimation error covariance ellipse for point F');

saveas(gcf, 'ellipse_F.eps', 'epsc');


