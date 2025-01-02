clc;clearvars

global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso N M_new M_new_2 theta_RIS status hncs gn h_SISO g_nkl
siso = 1; % 1 ise siso var, 0 ise yok
env = 1;
Frequency=28;
Noise = db2pow(-130);      % -100 dBm

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];
RIS_xyz=[40,50,1];

A = zeros(2*N,1); %inequality constraints
A(1:2:2*N) = -1;
A = diag(A);      % such that -z1<= 0, -z3 <= 0 and so on.
b = zeros(2*N,1); % it makes real part of every variable greater than 0
rng('shuffle');

f = ones(1,2*N);

options = optimoptions('fmincon');
options.Display = 'iter-detailed';
options.MaxFunctionEvaluations=Inf;

% H=zeros(1,trial);
% status = "ideal";
% for ii = 1:trial
%     tic
%     fprintf("interation no:" + ii);    
%     [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
%     H(ii) = channel_calculate(f);
%     toc
% end
% H_ort = -mean(H);

%% HORT FINDER
H_ort = 2.9690e-12;
status = "ideal";
H_calc = 0;
while(abs(H_ort*1e12-H_calc*1e12) > 0.01)
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    H_calc = -channel_calculate(f);
    disp(abs(H_ort*1e12-H_calc*1e12));
end

% hold = zeros(1,100);
% for kk = 1:100
%     status = "ideal";
%     H_calc = 0;
%     while(abs(H_ort*1e12-H_calc*1e12) > 0.1)
%         [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
%         H_calc = -channel_calculate(f);
%         disp(abs(H_ort*1e12-H_calc*1e12));
%     end
%     status = "realistic";
%     z0=0.5*ones(1,2*N);
%     options.MaxIterations = 0;
%     [Z_mincon, fval] = fmincon(@channel_calculate,z0,A,b,[],[],[],[],[],options);
%     hold(kk) = -fval;
%     disp(kk);
% end

%% BRING THE VARIABLES that satisfy the expected channel H
struct = load('./mat_files/it500feedw100.mat','hncs','gn','M_new','M_new_2','g_nkl','h_SISO','theta_RIS');
hncs = struct.hncs;
gn = struct.gn;
M_new = struct.M_new;
M_new_2 = struct.M_new_2;
g_nkl = struct.g_nkl;
theta_RIS = struct.theta_RIS;
h_SISO = struct.h_SISO;

%% OPTIMIZATION
status = "realistic";
H_ort = 2.9690e-12;
z0=0.5*ones(1,2*N);
options.MaxIterations = 30;
[Z_mincon, fval] = fmincon(@channel_calculate,z0,A,b,[],[],[],[],[],options);
save('./mat_files/SISOit30.mat')

z0=Z_mincon;
options.MaxIterations = 100;
[Z_mincon, fval] = fmincon(@channel_calculate,z0,A,b,[],[],[],[],[],options);
save('./mat_files/SISOit100feedw30.mat')

% z0=Z_mincon;
% options.MaxIterations = 500;
% [Z_mincon, fval] = fmincon(@channel_calculate,z0,A,b,[],[],[],[],[],options);
% save('./mat_files/SISOit500feedw100.mat')

fmin = -fval*10^(-16);
PtdBm=0:5:30;  % Tx power (dBm) - max: 30 dBm = 0 dBW
Pt=10.^((PtdBm-30)/10);  % in Watt
R_fmin = log2(1+Pt*fmin/Noise);
R_ideal = log2(1+Pt*H_ort/Noise);

plot(0:5:30,R_fmin,'r-o','MarkerFaceColor','r');
title('Achievable Rate vs. Transmit Power ')
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
grid on;
hold on;
plot(0:5:30,R_ideal,'b--o','MarkerFaceColor','b');
legend('Realistic Response','Ideal Response');
hold off;
