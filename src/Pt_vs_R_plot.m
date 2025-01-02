clc;clearvars;
global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso
tic  
siso = 1; % 1 ise siso var, 0 ise yok
env = 1; 
Frequency=28; 
pdfpoints = 0;
Num = 1000;
rng('shuffle');
Noise = db2pow(-130);      % -100 dBm

%% zris = 1, n = 256

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];  
RIS_xyz=[40,50,1];


Channel=zeros(1,Num);
SNR=Channel;
Capacity = zeros(1,11);

count=1;
for PtdBm=0:3:30  % Tx power (dBm) - max: 30 dBm = 0 dBW
    Pt=10^((PtdBm-30)/10);  % in Watts
    
    for repeat=1:Num
        Channel(repeat)= new_cha2();
        SNR(repeat)=Pt*abs(Channel(repeat))^2 / Noise;
    end

    Capacity(count)=mean(log2(1+SNR)); %achievable rate
    count=count+1;
    fprintf("zris1 256 PtdBm ="+PtdBm+"\n");
end

save("us_zris1_256.mat")

plot(0:3:30,Capacity,'r-o','MarkerFaceColor','r');
title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + Num + " N = " + N)
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
ylim([0 15])
xlim([0 31])
grid on;
hold on;

%% zris = 2, n = 256

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];  
RIS_xyz=[40,50,2];


Channel=zeros(1,Num);
SNR=Channel;
Capacity = zeros(1,11);

count=1;
for PtdBm=0:3:30  % Tx power (dBm) - max: 30 dBm = 0 dBW
    Pt=10^((PtdBm-30)/10);  % in Watts
    
    for repeat=1:Num
        Channel(repeat)= new_cha2();
        SNR(repeat)=Pt*abs(Channel(repeat))^2 / Noise;
    end

    Capacity(count)=mean(log2(1+SNR)); %achievable rate
    count=count+1;
    fprintf("zris2 256 PtdBm ="+PtdBm+"\n");
end

save("us_zris2_256.mat")

plot(0:3:30,Capacity,'r-^','MarkerFaceColor','r');
title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + Num + " N = " + N)
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
ylim([0 15])
xlim([0 31])
grid on;
hold on;

%% zris = 1, n = 1024

Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];  
RIS_xyz=[40,50,1];

N = 1024;

Channel=zeros(1,Num);
SNR=Channel;
Capacity = zeros(1,11);

count=1;
for PtdBm=0:3:30  % Tx power (dBm) - max: 30 dBm = 0 dBW
    Pt=10^((PtdBm-30)/10);  % in Watts
    
    for repeat=1:Num
        Channel(repeat)= new_cha2();
        SNR(repeat)=Pt*abs(Channel(repeat))^2 / Noise;
    end

    Capacity(count)=mean(log2(1+SNR)); %achievable rate
    count=count+1;
    fprintf("zris1 1024 PtdBm ="+PtdBm+"\n");
end

save("us_zris1_1024.mat")

plot(0:3:30,Capacity,'b-o','MarkerFaceColor','b');
title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + Num + " N = " + N)
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
ylim([0 15])
xlim([0 31])
grid on;
hold on;

%% zris = 2, n = 1024

Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];  
RIS_xyz=[40,50,2];

N = 1024;

Channel=zeros(1,Num);
SNR=Channel;
Capacity = zeros(1,11);

count=1;
for PtdBm=0:3:30  % Tx power (dBm) - max: 30 dBm = 0 dBW
    Pt=10^((PtdBm-30)/10);  % in Watts
    
    for repeat=1:Num
        Channel(repeat)= new_cha2();
        SNR(repeat)=Pt*abs(Channel(repeat))^2 / Noise;
    end

    Capacity(count)=mean(log2(1+SNR)); %achievable rate
    count=count+1;
    fprintf("zris2 1024 PtdBm ="+PtdBm+"\n");
end

save("us_zris2_1024.mat")

plot(0:3:30,Capacity,'b-^','MarkerFaceColor','b');
title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + Num + " N = " + N)
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
ylim([0 15])
xlim([0 31])
grid on;
hold on;









