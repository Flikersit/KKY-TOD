clc;
clear all;
close all;
%% === init===

T  = 1;
F  = [1 T; 0 1];           
H  = [1 0];                
Q  = (1/40) * [T^3/3  T^2/2;  
               T^2/2  T   ];
R  = 10;                    
Px0 = [10 0; 0 1];
x0m = [10; 1];

%% ===Kalmans filter===

syms x_00 [2 1] real
syms P_00 [2 2] real

syms x_10 [2 1] real
syms P_10 [2 2] real

syms K_1 [2 1] real
syms z real

%
x_10r = F*x_00;
P_10r = F*P_00*F' + Q;
% kalmans gain
K_1r = P_10*H'*inv(H*P_10*H' + R);
%
x_11 = x_10 + K_1*(z - H*x_10);
P_11 = (eye(2) - K_1*H)*P_10*(eye(2) - K_1*H)' + K_1*R*K_1';
%
inov_kalmn = z - H*x_10;


%% ===infinity kalmans filter===
syms x_10 [2 1] real
syms P_10 [2 2] real

syms z real

I = eye(size(F));

[Pp, pom, ~] = idare(F', H', Q, R, [], I);
Kinf = Pp*H'*inv(H*Pp*H' + R);

x_11_inf = x_10 + Kinf*(z - H*x_10);
P_11_inf = (eye(2) - Kinf*H)*P_10*(eye(2) - Kinf*H)' + Kinf*R*Kinf';


%% == State reconstructor (deterministic observer)==
syms x_00 [2 1] real
syms z_1 real

K = [1;
     1];

rekonstruktor = (eye(2) - K*H)*F*x_00 + K*z_1;  

%% ===SIMULATION===
Ntraj  = 1000;
Ksteps = 26;

P_filt_theory = zeros(2, 2, 26);
P_filt_inf_theory = zeros(2, 2, 26);

inovace_kalmn = zeros(Ntraj, Ksteps);
inovace_kalmn_inf = zeros(Ntraj, Ksteps);

x_true_store = zeros(Ntraj,2,Ksteps);

z_store = zeros(Ntraj, Ksteps);

x_pred_store = zeros(Ntraj,2,Ksteps);     % x_{k|k-1} Kalmn(classic)
x_filt_store = zeros(Ntraj,2,Ksteps);     % x_{k|k} Kalmn(classic)

x_pred_inf_store = zeros(Ntraj,2,Ksteps);     % x_{k|k-1} Kalmn(inf)
x_filt_inf_store = zeros(Ntraj,2,Ksteps);     % x_{k|k} Kalmn(inf)

x_rekonst_store = zeros(Ntraj,2,Ksteps);

for n = 1:Ntraj
    
    x_hat_rekonst = x0m;
    x_hat_inf = x0m;
    x_hat = x0m;
    
    x_true = mvnrnd(x0m, Px0)';
                  
    P     = Px0;                  
    P_inf = Px0;
    
    % --- inicializace  ---
    x_pred  = x_hat;   % x_{0|-1}
    P_pred  = P;       % P_{0|-1}
    x_pred_inf = x_hat_inf;
    P_pred_inf = P_inf;
    
    for k = 1:Ksteps
        
        w = mvnrnd([0;0], Q)';      
        v = sqrt(R)*randn(1);
        
        x_true = F*x_true + w;      % x_k
        z_k    = H*x_true + v;      % z_k
        
        z_store(n, k) = z_k;
        x_true_store(n,:,k) = x_true';
        
        %%%%%%%%%%%%%%%% FILTR %%%%%%%%%%%%%%%%%%
        
        % === Classic KF ===
        k1 = double(subs(K_1r, P_10, P_pred));
        
        tmp  = subs(inov_kalmn, x_10, x_pred);
        inov = double(subs(tmp, z, z_k));
        inovace_kalmn(n,k) = inov;
        
        tmp = subs(x_11, x_10, x_pred);
        tmp = subs(tmp, K_1, k1);
        x_filt = double(subs(tmp, z, z_k));
        
        tmpP = subs(P_11, P_10, P_pred);
        tmpP = subs(tmpP, K_1, k1);
        P_filt = double(tmpP);
        
        % === Rekonstruktor ===
        tmpr = subs(rekonstruktor, x_00, x_hat_rekonst);
        x_rek = double(subs(tmpr, z_1, z_k));
        
        
        % === Infinite KF ===
        tmp  = subs(inov_kalmn, x_10, x_pred_inf);
        inov = double(subs(tmp, z, z_k));
        inovace_kalmn_inf(n,k) = inov;
        
        tmp = subs(x_11_inf, x_10, x_pred_inf);
        x_filt_inf = double(subs(tmp, z, z_k));
        
        tmpP = subs(P_11_inf, P_10, P_pred_inf);
        P_filt_inf  = double(tmpP);
        
        
        %%%%%%%%%%%%%%%% PREDICTION %%%%%%%%%%%%%
        
        
        x_pred  = double(subs(x_10r, x_00, x_filt));
        P_pred  = double(subs(P_10r, P_00, P_filt));
        
        x_pred_inf  = double(subs(x_10r, x_00, x_filt_inf));
        P_pred_inf  = double(subs(P_10r, P_00, P_filt_inf));
        
        
        %%%%%%%%%%%% store %%%%%%%%%%%%
        
        x_pred_store(n,:,k) = x_pred';
        x_filt_store(n,:,k) = x_filt';
        
        x_pred_inf_store(n,:,k) = x_pred_inf';
        x_filt_inf_store(n,:,k) = x_filt_inf';
        
        x_rekonst_store(n,:,k) = x_rek';
        
        % update estimates
        x_hat_rekonst = x_rek;
        
        x_hat = x_filt;
        P     = P_filt;
        
        x_hat_inf = x_filt_inf;
        P_inf = P_filt_inf;
        
        if n == 1
           P_filt_theory(:,:,k)     = P_filt; 
           P_filt_inf_theory(:,:,k) = P_filt_inf;
        end
        
    end
end


%%
save('simulace.mat')

%%
load('simulace.mat')


%% ==== VISUALIZATION OF ONE TRAJECTORY ====

traj = 1;
t = 1:Ksteps;

figure;

subplot(2,1,1)
plot(t, squeeze(x_true_store(traj,1,:)), 'k-', 'LineWidth', 1.0); hold on;
plot(t, squeeze(x_true_store(traj,1,:)), 'kx', 'MarkerSize', 6);

plot(t, squeeze(x_pred_store(traj,1,:)), 'r--', 'LineWidth', 1.0);
plot(t, squeeze(x_pred_store(traj,1,:)), 'rx', 'MarkerSize', 6);

plot(t, squeeze(x_filt_store(traj,1,:)), 'r-', 'LineWidth', 1.0);
plot(t, squeeze(x_filt_store(traj,1,:)), 'rx', 'MarkerSize', 6);

plot(t, squeeze(x_pred_inf_store(traj,1,:)), 'b--', 'LineWidth', 1.0);
plot(t, squeeze(x_pred_inf_store(traj,1,:)), 'bx', 'MarkerSize', 6);

plot(t, squeeze(x_filt_inf_store(traj,1,:)), 'b-', 'LineWidth', 1.0);
plot(t, squeeze(x_filt_inf_store(traj,1,:)), 'bx', 'MarkerSize', 6);

plot(t, squeeze(x_rekonst_store(traj,1,:)), 'g-', 'LineWidth', 1.0);
plot(t, squeeze(x_rekonst_store(traj,1,:)), 'gx', 'MarkerSize', 6);

plot(t, z_store(1, :), 'cx', 'MarkerSize', 6);

grid on;
xlabel('k');
ylabel('x_1');
title('Position x_1');
legend('True state','True state (points)', ...
       'Prediction','Prediction (points)', ...
       'Filtering','Filtering (points)', ...
       'Prediction inf','Prediction inf (points)', ...
       'Filtering inf','Filtering inf (points)', ...
       'Observer', 'Observer(points)', ...
       'Measurements');

subplot(2,1,2)
plot(t, squeeze(x_true_store(traj,2,:)), 'k-', 'LineWidth', 1.0); hold on;
plot(t, squeeze(x_true_store(traj,2,:)), 'kx', 'MarkerSize', 6);

plot(t, squeeze(x_pred_store(traj,2,:)), 'r--', 'LineWidth', 1.0);
plot(t, squeeze(x_pred_store(traj,2,:)), 'rx', 'MarkerSize', 6);

plot(t, squeeze(x_filt_store(traj,2,:)), 'r-', 'LineWidth', 1.0);
plot(t, squeeze(x_filt_store(traj,2,:)), 'rx', 'MarkerSize', 6);

plot(t, squeeze(x_pred_inf_store(traj,2,:)), 'b--', 'LineWidth', 1.0);
plot(t, squeeze(x_pred_inf_store(traj,2,:)), 'bx', 'MarkerSize', 6);

plot(t, squeeze(x_filt_inf_store(traj,2,:)), 'b-', 'LineWidth', 1.0);
plot(t, squeeze(x_filt_inf_store(traj,2,:)), 'bx', 'MarkerSize', 6);

plot(t, squeeze(x_rekonst_store(traj,2,:)), 'g-', 'LineWidth', 1.0);
plot(t, squeeze(x_rekonst_store(traj,2,:)), 'gx', 'MarkerSize', 6);

plot(t, z_store(1, :), 'cx', 'MarkerSize', 6);

grid on;
xlabel('k');
ylabel('x_2');
title('Velocity x_2');
legend('True state','True state (points)', ...
       'Prediction','Prediction (points)', ...
       'Filtering','Filtering (points)', ...
       'Prediction inf','Prediction (points) inf', ...
       'Filtering inf','Filtering (points) inf', ...
       'Observer', 'Observer(points)', ...
       'Measurements');
saveas(gcf, 'first_traektory.eps', 'epsc');   
%% errors

err_kalmn = x_true_store - x_filt_store;
err_kalmn_inf = x_true_store - x_filt_inf_store;
err_rekonst = x_true_store - x_rekonst_store;


traj = 1;
t = 1:Ksteps;

figure;

subplot(2,1,1)
plot(t, squeeze(err_kalmn(traj,1,:)), 'r-', 'LineWidth', 1.0); hold on;
plot(t, squeeze(err_kalmn_inf(traj,1,:)), 'g-', 'LineWidth', 1.0);
plot(t, squeeze(err_rekonst(traj,1,:)), 'b-', 'LineWidth', 1.0); 

grid on;
xlabel('k');
ylabel('x_1');
title('Position error x_1');
legend('Error kalmn', 'Error inf kalman', 'Error observer');

subplot(2,1,2)
plot(t, squeeze(err_kalmn(traj,2,:)), 'r-', 'LineWidth', 1.0); hold on;
plot(t, squeeze(err_kalmn_inf(traj,2,:)), 'g-', 'LineWidth', 1.0);
plot(t, squeeze(err_rekonst(traj,2,:)), 'b-', 'LineWidth', 1.0);

grid on;
xlabel('k');
ylabel('x_2');
title('Velocity error x_2');
legend('Error kalmn', 'Error inf kalman', 'Error observer');


saveas(gcf, 'errors.eps', 'epsc'); 

%% === INNOVATION SEQUENCE PLOTS ===
figure;
plot(1:Ksteps, inovace_kalmn(1,:), 'LineWidth', 2);
hold on;
plot(1:Ksteps, inovace_kalmn_inf(1,:), 'LineWidth', 2);
grid on;
legend('Classic filter', 'Inf gain filter');
title('Comparison of innovation sequences');
xlabel('k');
ylabel('innovation');

saveas(gcf, 'innovations.eps', 'epsc'); 

%% === MEAN INNOVATION PLOTS ===
mu_kalmn = mean(inovace_kalmn, 1);
mu_kalmn_inf = mean(inovace_kalmn_inf, 1);

figure;
plot(1:Ksteps, mu_kalmn, 'bx', 'MarkerSize', 6);
hold on;
plot(1:Ksteps, mu_kalmn_inf, 'rx', 'MarkerSize', 6);
grid on;
legend('Classic filter', 'Inf gain filter');
title('Comparison of innovation sequences');
xlabel('k');
ylabel('innovation');

saveas(gcf, 'mean_innovations.eps', 'epsc'); 



%% === VARIANCES TABLE: C(k,k) for k = 0,1,5,25 ===
k_assignment = [0, 1, 5, 25];
k_indices = [1, 2, 6, 26]; % MATLAB indices
k_list = [1, 2, 5];
N = length(k_indices);
CovTab_kalmn = zeros(N, N);
CovTab_kalmn_inf = zeros(N, N);

for i = 1:N
    for j = 1:N
        Cmat = cov(inovace_kalmn(:, k_indices(i)), inovace_kalmn(:, k_indices(j)));
        CovTab_kalmn(i,j) = Cmat(1,2); % covariance between times
        Cmat = cov(inovace_kalmn_inf(:, k_indices(i)), inovace_kalmn_inf(:, k_indices(j)));
        CovTab_kalmn_inf(i,j) = Cmat(1,2);
    end
end

fprintf("\n===== Covariance Matrix of Innovations (Standart kalmn) =====\n");
fprintf("         k=0        k=1        k=5       k=25\n");
fprintf("--------------------------------------------------\n");
for i = 1:N
    fprintf("k=%-2d | ", k_assignment(i));
    fprintf("%8.4f  ", CovTab_kalmn(i,:));
    fprintf("\n");
end
fprintf("--------------------------------------------------\n");



fprintf("\n===== Covariance Matrix of Innovations (Inf kalmn) =====\n");
fprintf("         k=0        k=1        k=5       k=25\n");
fprintf("--------------------------------------------------\n");
for i = 1:N
    fprintf("k=%-2d | ", k_assignment(i));
    fprintf("%8.4f  ", CovTab_kalmn_inf(i,:));
    fprintf("\n");
end
fprintf("--------------------------------------------------\n");



%% === COVARIANCE BETWEEN TIMES C(k, k+l) ===
Ckl = cell(length(k_list),1);

fprintf('\nEstimated Covariance C(k, k+l)\n');
fprintf('----------------------------------------\n');
fprintf('  k    l    C(k,k+l)\n');
fprintf('----------------------------------------\n');

for idx = 1:length(k_list)
    k0 = k_list(idx);
    max_l = Ksteps - k0;
    C_row = zeros(1, max_l+1);

    for l = 0:max_l
        Cmat = cov(inovace_kalmn(:,k0), inovace_kalmn(:,k0+l));
        C_row(l+1) = Cmat(1,2);
        fprintf('%3d   %2d   %.5f\n', k_assignment(idx), l, C_row(l+1));
    end

    fprintf('----------------------------------------\n');
    Ckl{idx} = C_row;
end


%% === PLOTS: time-correlations C(k,k+l) ===
figure;
for idx = 1:length(k_list)
    C_row = Ckl{idx};
    l_values = 0:(length(C_row)-1);

    subplot(3,1,idx);
    stem(l_values, C_row, 'LineWidth', 1.8, 'Marker','o');
    grid on;
    yline(0,'k--');
    xlabel('l');
    ylabel('C(k,k+l)');
    title(['Covariance sequence for k = ', num2str(k_assignment(idx))]);
end

%% === COVARIANCE BETWEEN TIMES C(k, k+l) (INF)===
Ckl = cell(length(k_list),1);

fprintf('\nEstimated Covariance C(k, k+l)\n');
fprintf('----------------------------------------\n');
fprintf('  k    l    C(k,k+l)\n');
fprintf('----------------------------------------\n');

for idx = 1:length(k_list)
    k0 = k_list(idx);
    max_l = Ksteps - k0;
    C_row = zeros(1, max_l+1);

    for l = 0:max_l
        Cmat = cov(inovace_kalmn_inf(:,k0), inovace_kalmn_inf(:,k0+l));
        C_row(l+1) = Cmat(1,2);
        fprintf('%3d   %2d   %.5f\n', k_assignment(idx), l, C_row(l+1));
    end

    fprintf('----------------------------------------\n');
    Ckl{idx} = C_row;
end


%% === PLOTS: time-correlations C(k,k+l) (INF)===
figure;
for idx = 1:length(k_list)
    C_row = Ckl{idx};
    l_values = 0:(length(C_row)-1);

    subplot(3,1,idx);
    stem(l_values, C_row, 'LineWidth', 1.8, 'Marker','o');
    grid on;
    yline(0,'k--');
    xlabel('l');
    ylabel('C(k,k+l)');
    title(['Covariance sequence for k = ', num2str(k_assignment(idx))]);
end

%% === MSE ===
mse_kalmn_x1     = squeeze(mean(err_kalmn(:,1,:).^2,     1));   % [1 x Ksteps]
mse_kalmn_inf_x1 = squeeze(mean(err_kalmn_inf(:,1,:).^2, 1));
mse_rek_x1       = squeeze(mean(err_rekonst(:,1,:).^2,   1));


mse_kalmn_x2     = squeeze(mean(err_kalmn(:,2,:).^2,     1));
mse_kalmn_inf_x2 = squeeze(mean(err_kalmn_inf(:,2,:).^2, 1));
mse_rek_x2       = squeeze(mean(err_rekonst(:,2,:).^2,   1));

figure;

subplot(2,1,1);
hold on; grid on;

% МSE
plot(mse_kalmn_x1, 'b', 'LineWidth', 1.6);
plot(mse_kalmn_inf_x1, 'r', 'LineWidth', 1.6);
plot(mse_rek_x1, 'g', 'LineWidth', 1.6);

plot(squeeze(P_filt_theory(1,1,:)), 'b--', 'LineWidth', 1.2);
plot(squeeze(P_filt_inf_theory(1,1,:)), 'r--', 'LineWidth', 1.2);

title('Position x_1 MSE vs. theoretical covariance');
ylabel('MSE');
legend('KF','KF ∞-gain','Observer','KF theory','KF ∞ theory', ...
       'Location','best');

subplot(2,1,2);
hold on; grid on;

% МSE
plot(mse_kalmn_x2, 'b', 'LineWidth', 1.6);
plot(mse_kalmn_inf_x2, 'r', 'LineWidth', 1.6);
plot(mse_rek_x2, 'g', 'LineWidth', 1.6);

plot(squeeze(P_filt_theory(2,2,:)), 'b--', 'LineWidth', 1.2);
plot(squeeze(P_filt_inf_theory(2,2,:)), 'r--', 'LineWidth', 1.2);

title('Velocity x_2 MSE vs. theoretical covariance');
xlabel('k'); ylabel('MSE');
legend('KF','KF ∞-gain','Observer','KF theory','KF ∞ theory', ...
       'Location','best');
   
saveas(gcf, 'mse.eps', 'epsc'); 
