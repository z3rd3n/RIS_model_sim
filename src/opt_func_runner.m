clc;clearvars;
global env Tx_xyz RIS_xyz Rx_xyz N Frequency siso P
siso = 1; % 1 ise siso var, 0 ise yok
env = 1; 
Frequency=28; 
Noise = db2pow(-130);      % -100 dBm

N = 256;
Tx_xyz=[0,25,2];
Rx_xyz=[38,48,1];  
RIS_xyz=[40,50,1];

% A = zeros(2*N,1); %inequality constraints
% A(1:2:2*N) = -1;
% A = diag(A);      % such that -z1<= 0, -z3 <= 0 and so on.
% b = zeros(2*N,1); % it makes real part of every variable greater than 0

lb = zeros(1,2*N);
for dd = 1:length(lb)/2
    lb(2*dd) = -inf;
end

PtdBm=0:5:30;  % Tx power (dBm) - max: 30 dBm = 0 dBW
Pt=10.^((PtdBm-30)/10);  % in Watt

Capacity = zeros(1,length(Pt));

for pp = 1:length(Pt)
P = Pt(pp); 
[z, R] = particleswarm(@opt_func,2*N,lb);
Capacity(pp) = R;
end

plot(0:5:30,Capacity,'r-o','MarkerFaceColor','r');
title('Achievable Rate vs. Transmit Power ')
ylabel('Achievable Rate [bit/sec/Hz]')
xlabel('Transmit Power [dBm]')
grid on;
hold on;

