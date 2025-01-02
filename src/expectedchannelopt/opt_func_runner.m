clc;clearvars;
global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso 

siso = 0; % 1 ise siso var, 0 ise yok
env = 1;
Frequency=28;

Noise = db2pow(-130);      % -100 dBm
rng('shuffle');
N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];
RIS_xyz=[40,50,1];

A = zeros(2*N,1); %inequality constraints
A(1:2:2*N) = -1;
A = diag(A);      % such that -z1<= 0, -z3 <= 0 and so on.
b = zeros(2*N,1); % it makes real part of every variable greater than 0

z0=0.5*ones(1,2*N);

options = optimoptions('fmincon');
options.Display = 'iter-detailed';
options.BarrierParamUpdate = 'predictor-corrector';
options.MaxFunctionEvaluations=Inf;


PtdBm=0:5:30;  % Tx power (dBm) - max: 30 dBm = 0 dBW
Pt=10.^((PtdBm-30)/10);  % in Watt
a=zeros(1,10);
parfor mm = 1:10
    a(mm) = opt_func(ones(1,512));
end
%[Z_mincon, fval] = fmincon(@opt_func,z0,A,b,[],[],[],[],[],options);
fval = fval*1e-16;
SNR=Pt*-fval / Noise;
R=log2(1+SNR); %achievable rate

plot(0:5:30,R,'r-o','MarkerFaceColor','r');
title('Achievable Rate vs. Transmit Power ')
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
grid on;
hold on;

