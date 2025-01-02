clc;clear

global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso N M_new M_new_2 theta_RIS status hncs gn h_SISO g_nkl
siso = 1; % 1 ise siso var, 0 ise yok
env = 1;
Frequency=28;
Noise = db2pow(-130);      % -100 dBm

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];
RIS_xyz=[40,50,1];


rng('shuffle');
f = ones(1,2*N);

trial = 1;
hold_SO = zeros(1,trial);
hold_PS = zeros(1,trial);
hold_GA = zeros(1,trial);
hold_ideal = zeros(1,trial);


lb_SO = 1e-12*ones(1,N);
ub = 10e-12*ones(1,N);


options = optimoptions('surrogateopt');
% options.MinSurrogatePoints = max(20,6*N);
% options.MinSampleDistance = 2*1e-2;
% options.MaxFunctionEvaluations = 2250;
% options.MaxTime = 500;


for ii = 1:trial
    tic
    fprintf("interation no:" + ii);
    status = "realistic";
    
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    
    [F_SO, C_SO] = surrogateopt(@channel_calculateVAR,lb_SO,ub,options);
    hold_SO(ii) = -C_SO;
    
    [F_PS, C_PS] = particleswarm(@channel_calculateVAR,N);
    hold_PS(ii) = -C_PS;
    
    [F_GA, C_GA] = ga(@channel_calculateVAR,N);
    hold_GA(ii) = -C_GA;
    
    status = "ideal";
    hold_ideal(ii) = -channel_calculateVAR(f);
    
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

save("opt_VAR.mat");

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



