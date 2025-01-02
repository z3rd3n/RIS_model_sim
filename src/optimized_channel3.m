clc;clearvars

global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso N M_new M_new_2 theta_RIS status hncs gn h_SISO g_nkl
siso = 0; % 1 ise siso var, 0 ise yok
env = 1;
Frequency=28;
Noise = db2pow(-130);      % -100 dBm

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];
RIS_xyz=[40,50,1];


rng('shuffle');
f = ones(1,2*N);
lb_SO = zeros(1,2*N);
lb_PS = zeros(1,2*N);
ub = 1e9*ones(1,2*N);
for dd = 1:length(lb_SO)/2
    lb_SO(2*dd) = -1e9;
    lb_PS(2*dd) = -Inf;
end

A = zeros(2*N,1); %inequality constraints
A(1:2:2*N) = -1;
A = diag(A);      % such that -z1<= 0, -z3 <= 0 and so on.
b = zeros(2*N,1); % it makes real part of every variable greater than 0

trial = 100;
hold_SO = zeros(1,trial);
hold_PS = zeros(1,trial);
hold_GA = zeros(1,trial);
hold_ideal = zeros(1,trial);

options = optimoptions('surrogateopt');
options.MinSurrogatePoints = max(20,6*N);
options.MinSampleDistance = 2*1e-2;
options.MaxFunctionEvaluations = 2250;
options.MaxTime = 500;


for ii = 1:trial
    tic
    fprintf("interation no:" + ii);
    status = "realistic";
    
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    
    [F_SO, z_SO] = surrogateopt(@channel_calculate,lb_SO,ub,options);
    hold_SO(ii) = -z_SO;
    
    [F_PS, z_PS] = particleswarm(@channel_calculate,2*N,lb_PS);
    hold_PS(ii) = -z_PS;
    
    [F_GA, z_GA] = ga(@channel_calculate,2*N,A,b);
    hold_GA(ii) = -z_GA;
    
    status = "ideal";
    hold_ideal(ii) = -channel_calculate(f);
    
    toc
end

PtdBm=0:5:30;  % Tx power (dBm) - max: 30 dBm = 0 dBW
Pt=10.^((PtdBm-30)/10);  % in Watt
R_SO = zeros(1,length(Pt));
R_PS = zeros(1,length(Pt));
R_GA = zeros(1,length(Pt));
R_ideal = zeros(1,length(Pt));

for kk = 1:length(Pt)
    SNR_SO=(Pt(kk)*hold_SO / Noise);
    R_SO(kk)= mean(log2(1+SNR_SO)); %achievable rate
    
    SNR_PS=(Pt(kk)*hold_PS / Noise);
    R_PS(kk)= mean(log2(1+SNR_PS)); %achievable rate
    
    SNR_GA=(Pt(kk)*hold_GA / Noise);
    R_GA(kk)= mean(log2(1+SNR_GA)); %achievable rate
    
    SNR_ideal = Pt(kk)*hold_ideal/Noise;
    R_ideal(kk) = mean(log2(1+SNR_ideal));
end

save("mat_files/3opt_N256_noSISO_empedans_zris1.mat");

plot(0:5:30,R_ideal,'k--o','MarkerFaceColor','k');
hold on;
plot(0:5:30,R_SO,'b-o','MarkerFaceColor','b');
hold on;
plot(0:5:30,R_PS,'m-o','MarkerFaceColor','m');
hold on;
plot(0:5:30,R_GA,'r-o','MarkerFaceColor','r');
hold on;
title('Achievable Rate vs. Transmit Power ')
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
grid on;

legend('Ideal Response','Realistic Response - SO', 'Realistic Response - PS', 'Realistic Response - GA');



