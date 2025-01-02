% SIMRISv1 kullanilarak elde edilmis plotlardir.
% kullanilan kodlar asagida verilmistir
% 1000 ortalama alinmis, 200 nokta icin cizilmistir, siso yoktur env 2 secilmistir
% BU KODLARI VAROLAN BIR KODUN UZERINE CIZDIRIN LUTFEN
load("XvsR\v1ERT_N_64_XvsR.mat");
plot(linspace(20,70,TOTAL),R,'r-s','MarkerEdgeColor','b','MarkerFaceColor','k','MarkerIndices',round(linspace(1,TOTAL,26)));
hold on;
load("XvsR\v1ERT_N_256_XvsR.mat");
plot(linspace(20,70,TOTAL),R,'r-d','MarkerEdgeColor','b','MarkerFaceColor','k','MarkerIndices',round(linspace(1,TOTAL,26)));
hold on;
load("XvsR\v1ERT_N_1024_XvsR.mat");
plot(linspace(20,70,TOTAL),R,'r-o','MarkerEdgeColor','b','MarkerFaceColor','k','MarkerIndices',round(linspace(1,TOTAL,26)));
hold on;
title("MAKALE ILE ERTUGRUL HOCA'NIN KODLARI ARASINDAKI FARK")
grid on;
hold on;
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('x_{RIS} [m]')
xlim([19 71]);

% %% N = 64 %% ERT PLOTS
% N = 64;
% RIS_x = linspace(20,70,TOTAL);
% saved2 = zeros(1,TOTAL);
% saved = zeros(1,TOTAL2);
% for mm = 1:TOTAL
%     RIS_xyz(1) = RIS_x(mm);
%     for ii = 1:TOTAL2
%         fprintf("ertN = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
%         saved(ii) = SimRISv1();
%     end
%     saved2(mm) = mean(saved);
% end
% Pt = 1;
% No = db2pow(-130);
% SNR = (((abs(saved2).^2)/No)).*Pt;
% R = log2(1+SNR);
% save("v1ERT_N_64_XvsR.mat");
% plot(linspace(20,70,TOTAL),R,'b-s','MarkerEdgeColor','b','MarkerFaceColor','#D95319','MarkerIndices',round(linspace(1,TOTAL,26)));
% title("X_{RIS}  vs. R ")
% grid on;
% hold on;
% 
% 
% %% N = 256 %% ERT PLOTS
% N = 256;
% RIS_x = linspace(20,70,TOTAL);
% saved2 = zeros(1,TOTAL);
% saved = zeros(1,TOTAL2);
% for mm = 1:TOTAL
%     RIS_xyz(1) = RIS_x(mm);
%     for ii = 1:TOTAL2
%         fprintf("ertN = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
%         saved(ii) = SimRISv1();
%     end
%     saved2(mm) = mean(saved);
% end
% Pt = 1;
% No = db2pow(-130);
% SNR = (((abs(saved2).^2)/No)).*Pt;
% R = log2(1+SNR);
% save("v1ERT_N_256_XvsR.mat");
% plot(linspace(20,70,TOTAL),R,'b-d','MarkerEdgeColor','b','MarkerFaceColor','y','MarkerIndices',round(linspace(1,TOTAL,26)));
% title("X_{RIS}  vs. R ")
% grid on;
% hold on;
% 
% 
% %% N = 1024 %% ERT PLOTS
% N = 1024;
% RIS_x = linspace(20,70,TOTAL);
% saved2 = zeros(1,TOTAL);
% saved = zeros(1,TOTAL2);
% for mm = 1:TOTAL
%     RIS_xyz(1) = RIS_x(mm);
%     for ii = 1:TOTAL2
%         fprintf("ertN = " + N + "; mm = " + mm + "; ii = "+ii+"\n");
%         saved(ii) = SimRISv1();
%     end
%     saved2(mm) = mean(saved);
% end
% Pt = 1;
% No = db2pow(-130);
% SNR = (((abs(saved2).^2)/No)).*Pt;
% R = log2(1+SNR);
% save("v1ERT_N_1024_XvsR.mat");
% plot(linspace(20,70,TOTAL),R,'b-o','MarkerEdgeColor','b','MarkerFaceColor','g','MarkerIndices',round(linspace(1,TOTAL,26)));
% title("X_{RIS}  vs. R ")
% grid on;
% hold on;