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
status = "ideal";

rng('shuffle');
f = zeros(1,2*N);
trial = 1e4;
hold_ideal = zeros(1,trial);

for ii = 1:trial
    tic
    fprintf("interation no:" + ii);
    [M_new,M_new_2,hncs,gn,g_nkl,h_SISO,theta_RIS] = new_cha3();
    hold_ideal(ii) = -channel_calculate(f);
    
    toc
end

histogram(sqrt(aaa)*1e5,128,'Normalization','pdf');

pd = histfit((sqrt(aaa)*1e5)',128,'birnbaumsaunders');

