load("XvsR\n64.mat");
n64 = Data001;
load("XvsR\n256.mat");
n256 = Data002;
load("XvsR\n1024.mat");
n1024 = Data003;

plot(n64(:,1),n64(:,2),'b-s','MarkerEdgeColor','b','MarkerFaceColor','#D95319');
grid on
hold on
plot(n256(:,1),n256(:,2),'b-d','MarkerEdgeColor','b','MarkerFaceColor','y');
grid on
hold on
plot(n1024(:,1),n1024(:,2),'b-o','MarkerEdgeColor','b','MarkerFaceColor','g');
grid on
hold on
xlim([19 71]);
