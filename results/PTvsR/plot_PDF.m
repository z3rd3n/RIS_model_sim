% bu grafikler positioning makalesindeki grafikleri ifade eder
load("PTvsR\no_RIS.mat");
plot(no_RIS(:,1),no_RIS(:,2),'k-');
hold on
load("PTvsR\256zris1.mat");
zris_256_1 = Data001;
plot(zris_256_1(:,1),zris_256_1(:,2),'r:o','MarkerFaceColor','r');
hold on
load("PTvsR\256zris2.mat");
zris_256_2 = Data002;
plot(zris_256_2(:,1),zris_256_2(:,2),'r:^','MarkerFaceColor','r');
hold on
load("PTvsR\1024zris1.mat");
zris_1024_1 = Data003;
plot(zris_1024_1(:,1),zris_1024_1(:,2),'b:o','MarkerFaceColor','b');
hold on
load("PTvsR\1024zris2.mat");
zris_1024_2 = Data004;
plot(zris_1024_2(:,1),zris_1024_2(:,2),'b:^','MarkerFaceColor','b');
hold on
ylim([0 15])
xlim([0 31])
grid on;
hold on;



