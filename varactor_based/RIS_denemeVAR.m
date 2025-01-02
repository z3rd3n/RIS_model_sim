function slack = RIS_denemeVAR(ANGLE,C_1)
%d is the thickness of the bulky medium.
%mu is the vacuum permability
%e_r is relative electric permittivity
%D is periodicity
%W is the gap between patches 
%NOTE THAT W and D may vary in different directions -> Wx, Wy, Dx, Dy

frequency = 28;
mu = 4*pi*1e-7;
e_r=4.4-1i*0.088;
d = 1.2e-3;
D = 5e-3;
W = 0.5e-3;
sigma = 58.7*1e6;
Rvar = 0.5;
Lvar = 0.5e-9;
tech = "Single Polarized";

e_0 = 8.854*1e-12;
eta_0 = 376.730; %Free Space Impedance
mu_0 = 4*pi*1e-7; %vacuum permeability 
freq = frequency*1e9;
w = 2*pi*freq; % Angular Frequency
k = w/(3*10^8); % Wave Number
k_z = k*sqrt(e_r - sin(ANGLE)^2);% normal propagation constant within the substrate
delta = 2/sqrt(sigma*w*mu_0); %skin depth
e_eff = (1+e_r)/2; %effective electric permittivity
k_eff = k*sqrt(e_eff); %effective wavenumber
ANGLE_t = asin(sqrt(1/e_r)*sin(ANGLE)); %transmission angle

Z_d_TE = (1i*w*mu/k_z)*tan(k_z*d);
Z_d_TM = (k_z/(w*e_0*e_r))*tan(k_z*d);

R_patch = ((D/(D-W))^2)*(1/(sigma*delta)); %patch resistance

Dx = D; % There parameters can be adjusted if patches are not square!!!
Dy = D;
Wx = W;
Wy = W;

Ax = (-2*Dx*e_eff*e_0/pi)*log(sin(pi*Wx*0.5/Dx));
Ay = (-2*Dy*e_eff*e_0/pi)*log(sin(pi*Wy*0.5/Dy));

C_patch_TE = Ay*(1-(k/k_eff)*0.5*sin(ANGLE_t)^2);
C_patch_TM = Ax;

%ground-patch capacitance
Cap_correction_x=2*e_0*Dx/pi.*log(1-exp(-4*pi*d./Dx));
Cap_correction_y=2*e_0*Dy/pi.*log(1-exp(-4*pi*d./Dx));

Z_patch_TE = R_patch + 1/(1i*w*(C_patch_TE-Cap_correction_y));
Z_patch_TM = R_patch + 1/(1i*w*(C_patch_TM-Cap_correction_x));

Z_pi_TE = Z_d_TE*Z_patch_TE;
Z_pi_TM = Z_d_TM*Z_patch_TM;
Z_sm_TE = Z_d_TE+Z_patch_TE;
Z_sm_TM = Z_d_TM+Z_patch_TM;

eta_TE = (Z_sm_TE/Z_pi_TE)*eta_0;
eta_TM = (Z_sm_TM/Z_pi_TM)*eta_0;

Z_1 = Rvar + 1i*w*Lvar - 1i/(C_1*w);

if tech == "Dual Polarized"
    slack = [(Z_1*(1-eta_TE)-eta_0)/(Z_1*(1+eta_TE)+eta_0) 0;0 (Z_2*(1-eta_TM)-eta_0)/(Z_2*(1+eta_TM)+eta_0)];
elseif tech == "Single Polarized"
    slack = [(Z_1*(1-eta_TE)-eta_0)/(Z_1*(1+eta_TE)+eta_0) 0;0 (Z_1*(1-eta_TM)-eta_0)/(Z_1*(1+eta_TM)+eta_0)];
else
    error('You should specify a RIS technology!!!');
end

end