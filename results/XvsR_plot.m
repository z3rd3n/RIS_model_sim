clc;clearvars;
global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso
  
siso = 0; % 1 ise siso var, 0 ise yok
env = 2; 
Frequency = 28; 

TOTAL = 50;
TOTAL2 = 1000;

% recommeneded pozitions
Tx_xyz=[0,25,20];
Rx_xyz=[50,50,1];  
RIS_xyz=[20,60,10];

%% N = 64
N = 64;
RIS_x = linspace(20,70,TOTAL);
saved2 = zeros(1,TOTAL);
saved = zeros(1,TOTAL2);
for mm = 1:TOTAL
    RIS_xyz(1) = RIS_x(mm);
    for ii = 1:TOTAL2
        fprintf("N = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
        saved(ii) = new_cha2();
    end
    saved2(mm) = mean(saved);
end
Pt = 1;
No = db2pow(-130);
SNR = (((abs(saved2).^2)/No)).*Pt;
R = log2(1+SNR);
save("N_64_XvsR_" + datestr(date) + ".mat");
plot(linspace(20,70,TOTAL),R,'k-^','MarkerEdgeColor','k','MarkerFaceColor','#D95319','MarkerIndices',round(linspace(1,TOTAL,26)));
title("X_{RIS}  vs. R ")
grid on;
hold on;


%% N = 256
N = 256;
RIS_x = linspace(20,70,TOTAL);
saved2 = zeros(1,TOTAL);
saved = zeros(1,TOTAL2);
for mm = 1:TOTAL
    RIS_xyz(1) = RIS_x(mm);
    for ii = 1:TOTAL2
        fprintf("N = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
        saved(ii) = new_cha2();
    end
    saved2(mm) = mean(saved);
end
Pt = 1;
No = db2pow(-130);
SNR = (((abs(saved2).^2)/No)).*Pt;
R = log2(1+SNR);
save("N_256_XvsR_" + datestr(date) + ".mat");
plot(linspace(20,70,TOTAL),R,'k-v','MarkerEdgeColor','k','MarkerFaceColor','y','MarkerIndices',round(linspace(1,TOTAL,26)));
title("X_{RIS}  vs. R ")
grid on;
hold on;


%% N = 1024
N = 1024;
RIS_x = linspace(20,70,TOTAL);
saved2 = zeros(1,TOTAL);
saved = zeros(1,TOTAL2);
for mm = 1:TOTAL
    RIS_xyz(1) = RIS_x(mm);
    for ii = 1:TOTAL2
        fprintf("N = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
        saved(ii) = new_cha2();
    end
    saved2(mm) = mean(saved);
end
Pt = 1;
No = db2pow(-130);
SNR = (((abs(saved2).^2)/No)).*Pt;
R = log2(1+SNR);
save("N_1024_XvsR_" + datestr(date) + ".mat");
plot(linspace(20,70,TOTAL),R,'k-h','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerIndices',round(linspace(1,TOTAL,26)));
title("X_{RIS}  vs. R ")
grid on;
hold on;



