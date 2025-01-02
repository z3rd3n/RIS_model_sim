clc;clear

global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso N M_new M_new_2 theta_RIS status hncs gn h_SISO g_nkl
siso = 1; % 1 ise siso var, 0 ise yok
env = 1;
Frequency=28;
status = "realistic";
Noise = db2pow(-130);      % -100 dBm

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];
RIS_xyz=[40,50,1];


rng('shuffle');
f = ones(1,2*N);
lb = zeros(1,2*N);
ub = 1e9*ones(1,2*N);
for dd = 1:length(lb)/2
    lb(2*dd) = -1e9;
end

trial = 100;
hold_z = zeros(1,trial);
hold_ideal = zeros(1,trial);
options = optimoptions('surrogateopt');
options.MinSurrogatePoints = max(20,6*N);
options.MinSampleDistance = 2*1e-2;
options.MaxFunctionEvaluations = 2250;
options.MaxTime = 500;
for ii = 1:trial
    tic
    fprintf("iteration no:" + ii)
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    [F, z] = surrogateopt(@channel_calculate,lb,ub,options);
    hold_z(ii) = z;
    toc
end

z = -hold_z;

status = "ideal";

for jj = 1:trial
    fprintf("iteration no:" + jj)
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    hold_ideal(jj) = channel_calculate(f);
end

ideal_z = -hold_ideal;

PtdBm=0:5:30;  % Tx power (dBm) - max: 30 dBm = 0 dBW
Pt=10.^((PtdBm-30)/10);  % in Watt
R = zeros(1,length(Pt));
R_ideal = zeros(1,length(Pt));
for kk = 1:length(Pt)
    SNR=(Pt(kk)*z / Noise);
    R(kk)= mean(log2(1+SNR)); %achievable rate
    
    SNR_ideal = Pt(kk)*ideal_z/Noise;
    R_ideal(kk) = mean(log2(1+SNR_ideal));
end

plot(0:5:30,R,'c-o','MarkerFaceColor','c');
title('Achievable Rate vs. Transmit Power ')
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
grid on;
hold on;
plot(0:5:30,R_ideal,'g--o','MarkerFaceColor','g');
legend('Realistic Response','Ideal Response');

