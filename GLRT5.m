close all;
clear all;
%--------------------------------------------------------------------------
% Loading the dataset
load('Set1_hydrophone_data.mat');
load('Set2_hydrophone_data.mat');
data_1 = received_data_set1; 
data_2 = received_data_set2;
x_1 = mean(data_1,2); % Averaging the data to reduce noise
x_2 = mean(data_2,2); % Averaging the data to reduce noise
x_detect = x_2; % Change for dataset (x_1)
%--------------------------------------------------------------------------
% Given signal parameters
tr = 0.1;
fs = 800;
Ts = 1/fs;
Td = 1.5;
Trec = 4;
N = (Trec/Ts)+1;
t = linspace(0, N-1, N); 
M = ((Td + tr)/Ts)+1;
%--------------------------------------------------------------------------
% Range of Zeta Parameters 
zeta_min = 2;
zeta_max = 40;
n_zeta = zeta_min:0.1:zeta_max;
%--------------------------------------------------------------------------
% Range of delay values
n_min = 1;
n_max = N-M; 
n_o = n_min:1:n_max;
%--------------------------------------------------------------------------
% MLE values of zeta and delay
disp('Loop started');
for j = 1:length(n_zeta)
    for k = 1:length(n_o)
        c_zeta = zeros(1, N);
        s_zeta = zeros(1, N);
        H = zeros(N, 2);
            c_zeta(n_o(k):n_o(k)+M-1) = cos(2*pi*n_zeta(j)*log((tr+(0:1280)*Ts)/tr));
            s_zeta(n_o(k):n_o(k)+M-1) = sin(2*pi*n_zeta(j)*log((tr+(0:1280)*Ts)/tr));
        H = [c_zeta', s_zeta'];
        J(k,j) = x_detect'*H*(inv(H'*H))*(H')*x_detect;
    end
    disp(j);
end
disp('Loop ended');
figure(1),mesh(J),title('Maximization Function for \zeta MLE and n_o MLE'),xlabel('Delay n_o'),ylabel('\zeta index'),zlabel('Value of J(n)');
max_value = max(max(J));
[n_index, zeta_index] = find(J==max_value);
zeta_MLE = n_zeta(zeta_index);
n_MLE = n_o(n_index);
disp('zeta_MLE='),disp(zeta_MLE);
disp('n_o_MLE='),disp(n_MLE);
%--------------------------------------------------------------------------
% MLE values for Amplitude and Phase
c_zeta_MLE = zeros(1, N);
s_zeta_MLE = zeros(1, N);
c_zeta_MLE(n_MLE:(n_MLE+M-1)) = cos(2*pi*zeta_MLE*log((tr+(0:M-1)*Ts)/tr));
s_zeta_MLE(n_MLE:(n_MLE+M-1)) = sin(2*pi*zeta_MLE*log((tr+(0:M-1)*Ts)/tr));
H_MLE = [c_zeta_MLE' s_zeta_MLE'];
theta_MLE = (inv(H_MLE'*H_MLE))*(H_MLE')*x_detect;
A_MLE = sqrt((theta_MLE(1)^2) + (theta_MLE(2)^2));
P_MLE = atan(-theta_MLE(2)/theta_MLE(1));
disp('Amplitude_MLE='),disp(A_MLE);
disp('Phase_MLE='),disp(P_MLE);
%--------------------------------------------------------------------------
% MLE for Variances under H0 and H1
var_1 = (1/N)*(x_detect - (H*theta_MLE))'*(x_detect - (H*theta_MLE));
var_0 = (1/N)*(x_detect)'*(x_detect);
disp('Variance under H1='),disp(var_1);
disp('Variance under H0='),disp(var_0);
%--------------------------------------------------------------------------
% Decision statistic
zeta_test = zeta_min:0.1:zeta_max;
n_test = n_min:1:n_max;
disp('Loop Started');
for ii = 1:length(zeta_test)
    for jj = 1:length(n_test)
        y_test = zeros(N, 1);
        y_test(n_test(jj):(n_test(jj)+M-1)) = A_MLE*cos((2*pi*zeta_test(ii)*log(1+(0:M-1)*Ts/tr))+P_MLE);
        H1_test = (1/N)*(x_detect - y_test)'*(x_detect - y_test);
        H0_test = (1/N)*(x_detect')*(x_detect) ;
        t_x(jj, ii) = (H1_test/H0_test)^(-N/2);
    end
    disp(ii);
end
disp('Loop Ended');
figure(2),mesh(t_x),title('Test Statistic'),xlabel('Delay n_o'),ylabel('\zeta index'),zlabel('Test Value');
%--------------------------------------------------------------------------
% Calculating range and doppler scale, if target is present
target_range = (Ts*(n_MLE - 1))*(1500/2);
doppler_b = exp(P_MLE/(2*pi*zeta_MLE));