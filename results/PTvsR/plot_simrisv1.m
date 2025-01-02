
% SIMRISv1 kullanilarak elde edilmis plotlardir.
% kullanilan kodlar asagida verilmistir
% 10000 ortalama alinmis, siso vardir env 1 secilmistir
% BU KODLARI VAROLAN BIR KODUN UZERINE CIZDIRIN LUTFEN

load("PTvsR\v1ERT_N256zris1_25-Aug-2022.mat");
plot(pow2db(Pt)+30,R,'r--o','MarkerFaceColor','r','MarkerIndices',1:7:50);
load("PTvsR\v1ERT_N256zris2_25-Aug-2022.mat");
plot(pow2db(Pt)+30,R,'r--^','MarkerFaceColor','r','MarkerIndices',1:7:50);
load("PTvsR\v1ERT_N1024zris1_25-Aug-2022.mat");
plot(pow2db(Pt)+30,R,'b--o','MarkerFaceColor','b','MarkerIndices',1:7:50);
load("PTvsR\v1ERT_N1024zris2_25-Aug-2022.mat");
plot(pow2db(Pt)+30,R,'b--^','MarkerFaceColor','b','MarkerIndices',1:7:50);
ylim([0 15])
xlim([0 31])
grid on;
hold on;

% %% ERT zris = 1, n = 256
% 
% Tx_xyz=[0,25,2];
% Rx_xyz=[38,48,1];  
% RIS_xyz=[40,50,1];
% 
% N = 256;
% 
% vals = zeros(1,TOTAL);
% for ii = 1:TOTAL
%     fprintf("ERTzris1 256 ii ="+ii+"\n");
%     vals(ii) = SimRISv1();
% end
% saved = mean(vals);
% 
% 
% Pt = logspace(-3,0,50); % 1e-3 to 1
% No = db2pow(-130);      % -100 dBm
% 
% SNR = (((abs(saved)^2)/No)).*Pt;
% R = log2(1+SNR);
% save("v1ERT_N256zris1_" + datestr(date) + ".mat");
% plot(pow2db(Pt)+30,R,'k--o','MarkerFaceColor','k','MarkerIndices',1:7:50);
% title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + TOTAL + " N = " + N)
% ylabel('Achievable Rate [bit/sec/Hz]')
% xlabel('Transmit Power [dBm]')
% ylim([0 15])
% xlim([0 31])
% grid on;
% hold on;
% 
% %% ERT zris = 2, n = 256
% 
% Tx_xyz=[0,25,2];
% Rx_xyz=[38,48,1];  
% RIS_xyz=[40,50,2];
% 
% N = 256;
% 
% vals = zeros(1,TOTAL);
% for ii = 1:TOTAL
%     fprintf("zris2 256 ii ="+ii+"\n");
%     vals(ii) = SimRISv1();
% end
% saved = mean(vals);
% 
% 
% Pt = logspace(-3,0,50);
% No = db2pow(-130);
% 
% SNR = (((abs(saved)^2)/No)).*Pt;
% R = log2(1+SNR);
% 
% save("v1ERT_N256zris2_" + datestr(date) + ".mat");
% plot(pow2db(Pt)+30,R,'k--^','MarkerFaceColor','k','MarkerIndices',1:7:50);
% title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + TOTAL + " N = " + N)
% ylabel('Achievable Rate [bit/sec/Hz]')
% xlabel('Transmit Power [dBm]')
% %ylim([0 15])
% xlim([0 31])
% grid on;
% hold on;
% 
% %% ERT zris = 1, n = 1024
% 
% Tx_xyz=[0,25,2];
% Rx_xyz=[38,48,1];  
% RIS_xyz=[40,50,1];
% 
% N = 1024;
% 
% vals = zeros(1,TOTAL);
% for ii = 1:TOTAL
%     fprintf("zris1 1024 ii ="+ii+"\n");
%     vals(ii) = SimRISv1();
% end
% saved = mean(vals);
% 
% 
% Pt = logspace(-3,0,50);
% No = db2pow(-130);
% 
% SNR = (((abs(saved)^2)/No)).*Pt;
% R = log2(1+SNR);
% 
% save("v1ERT_N1024zris1_" + datestr(date) + ".mat");
% plot(pow2db(Pt)+30,R,'m--o','MarkerFaceColor','m','MarkerIndices',1:7:50);
% title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + TOTAL + " N = " + N)
% ylabel('Achievable Rate [bit/sec/Hz]')
% xlabel('Transmit Power [dBm]')
% %ylim([0 15])
% xlim([0 31])
% grid on;
% hold on;
% 
% %% ERT zris = 2, n = 1024
% 
% Tx_xyz=[0,25,2];
% Rx_xyz=[38,48,1];  
% RIS_xyz=[40,50,2];
% 
% N = 1024;
% 
% vals = zeros(1,TOTAL);
% for ii = 1:TOTAL
%     fprintf("zris2 1024 ii ="+ii+"\n");
%     vals(ii) = SimRISv1();
% end
% saved = mean(vals);
% 
% 
% Pt = logspace(-3,0,50);
% No = db2pow(-130);
% 
% SNR = (((abs(saved)^2)/No)).*Pt;
% R = log2(1+SNR);
% 
% save("v1ERT_N1024zris2_" + datestr(date) + ".mat");
% plot(pow2db(Pt)+30,R,'m--^','MarkerFaceColor','m','MarkerIndices',1:7:50);
% title('Achievable Rate vs. Transmit Power ' +" env = " + env + " ortalama sayisi = " + TOTAL + " N = " + N)
% ylabel('Achievable Rate [bit/sec/Hz]')
% xlabel('Transmit Power [dBm]')
% %ylim([0 15])
% xlim([0 31])
% grid on;
% hold on;
