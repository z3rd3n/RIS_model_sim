function H = new_cha2()
% Serden Sait Eranil
% Ramazan Umut Aktas
global k dis env N Tx_xyz RIS_xyz Rx_xyz Frequency siso

lambda=(3*10^8)/(Frequency*10^9);  % Wavelength
k=2*pi/lambda;                     % Wavenumber
dis=lambda/2;                      % RIS element spacing (can be modified to observe its effect on h and g)


x_Tx=Tx_xyz(1);y_Tx=Tx_xyz(2);z_Tx=Tx_xyz(3);
x_Rx=Rx_xyz(1);y_Rx=Rx_xyz(2);z_Rx=Rx_xyz(3);
x_RIS=RIS_xyz(1);y_RIS=RIS_xyz(2);z_RIS=RIS_xyz(3);


if env==1 % INDOORS
    % env 1 (InH - Office) - NLOS
    n_NLOS=3.19;          % Path Loss Exponent (Indoor Office NLOS)
    sigma_NLOS=8.29;      % Shadow Fading Term (dB) (Indoor Office NLOS)
    b_NLOS=0.06;          % Path Loss Parameter (Indoor Office NLOS)
    f0=24.2;              % Path Loss Parameter (GHz) (Indoor Office NLOS)
    xpr_mu_NLOS = 10;      % cross correlation mean values (DB)
    xpr_std_NLOS = 4;
    
    % env 1 (InH - Office) - LOS
    n_LOS=1.73;           % Path Loss Exponent (Indoor Office LOS)
    sigma_LOS=3.02;       % Shadow Fading Term (dB) (Indoor Office LOS)
    b_LOS=0;              % Path Loss Parameter (Indoor Office NLOS)
    xpr_mu_LOS = 11;
    xpr_std_LOS = 4;
    
else  % OUTDOORS
    % Enviroment 2 (UMi - Street Canyon) - NLOS
    n_NLOS=3.19;          % Path Loss Exponent (UMi Street Canyon NLOS)
    sigma_NLOS=8.2;       % Shadow Fading Term (dB) (UMi Street Canyon NLOS)
    b_NLOS=0;             % Path Loss Parameter (UMi Street Canyon NLOS)
    f0=24.2;              % Useless since b_NLOS=0
    xpr_mu_NLOS = 8;
    xpr_std_NLOS = 3;
    
    % Enviroment 2 (UMi - Street Canyon) - LOS
    n_LOS=1.98;          % Path Loss Exponent (UMi Street Canyon LOS)
    sigma_LOS=3.1;       % Shadow Fading Term (dB) (UMi Street Canyon LOS)
    b_LOS=0;             % Path Loss Parameter (UMi Street Canyon LOS)
    xpr_mu_LOS = 9;
    xpr_std_LOS =3;
end

% Parameter of Number of Clusters
if Frequency==28
    lambda_p=1.8;
elseif Frequency==73
    lambda_p=1.9;
end

d_T_RIS = norm(Tx_xyz-RIS_xyz);    % or sqrt(sum((RIS_xyz-Tx_xyz).^2))


% LOS Probability is Relatively Low for Indoors if d_T_RIS > 20
if env==1
    if z_RIS<z_Tx   % for ground level RIS
        % InH LOS Probability
        if d_T_RIS<= 1.2
            p_LOS=1;
        elseif 1.2<d_T_RIS && d_T_RIS<6.5
            p_LOS=exp(-(d_T_RIS-1.2)/4.7);
        else
            p_LOS=0.32*exp(-(d_T_RIS-6.5)/32.6);
        end
        
        I_LOS=randsrc(1,1,[1,0;p_LOS 1-p_LOS]);
        
    elseif z_RIS>=z_Tx % for an RIS mounted at a high place (100% LOS)
        I_LOS=1;
        
    end
    
elseif env==2
    % UMi LOS Probability
    p_LOS=min([20/d_T_RIS,1])*(1-exp(-d_T_RIS/39)) + exp(-d_T_RIS/39);
    
    I_LOS=randsrc(1,1,[1,0;p_LOS 1-p_LOS]);
    
end



if I_LOS==1 %LOS ACILARINI HESAPLA
    
    I_phi=sign(x_RIS-x_Tx);
    phi_RIS_LOS = I_phi* atand ( abs( x_RIS-x_Tx) / abs(y_RIS-y_Tx) );
    
    I_theta=sign(z_Tx-z_RIS);
    theta_RIS_LOS=I_theta * asind ( abs (z_RIS-z_Tx )/ d_T_RIS );
    
    
    % Link Attentuation (LOS)
    
    % Note: This is different than FSPL (Shadowing/Waveguiding effect included with n < 2)
    L_dB_LOS=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_RIS)- randn*sigma_LOS;
    hsiso_val = randn*sigma_LOS;
    L_LOS=10^(L_dB_LOS/10);
    
    Ic = 0;
    %%%%%%%%%%%%
else
    h_LOS=0;
    phi_RIS_LOS = 0;
    theta_RIS_LOS = 0;
end

for generate=1:100  % To ensure that at least one scatterer exist
    
    if I_LOS == 1
        C=max([1,poissrnd(lambda_p)])+1;  % Poisson distributed
        S=randi(30,1,C);
        S(1)=1;
    else
        C=max([1,poissrnd(lambda_p)]);  % Poisson distributed
        S=randi(30,1,C);
    end
    
    phi_Tx=[ ];
    theta_Tx=[ ];
    phi_av=zeros(1,C);
    theta_av=zeros(1,C);
    
    for counter=1:C
        phi_av(counter)  = rand*180-90;     % mean azimuth
        theta_av(counter)= rand*90-45;      % mean elevation  (needs update & most of the clusters are above ceiling)
        % for outdoors this might not be a big concern
        % cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
        phi_Tx   = [phi_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + phi_av(counter)];
        theta_Tx = [theta_Tx,log(rand(1,S(counter))./rand(1,S(counter)))*sqrt(25/2) + theta_av(counter)];
    end  % LAPLACIAN ANLGES TOWARD CSTH CLUSTERS
    % Cluster Distances
    a_c=1+rand(1,C)*(d_T_RIS-1);    % Cluster distances uniform [1,d_T_RIS] % can be modified later
    % Correction on Cluster Locations for Indoors
    
    if env ==1
        % Room dimensions (Indoor Hotspot)
        dim=[75,50,3.5];                  % x-y dimensions recommended by 5G Channel Model, height is assumed as 3.5 m
        % Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
        Coordinates=zeros(C,3);          % for Clusters
        Coordinates2=zeros(sum(S),3);    % for Scatterers
        for counter=1:C
            loop=1;
            Coordinates(counter,:)=[x_Tx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),... %scatterer positions determined
                y_Tx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                z_Tx + a_c(counter)*sind(theta_av(counter))] ;
            while Coordinates(counter,3)>dim(3) || Coordinates(counter,3)<0 ||  Coordinates(counter,2)>dim(2) ||  Coordinates(counter,2)<0  ||  Coordinates(counter,1)>dim(1) ||  Coordinates(counter,1)<0
                a_c(counter)=    0.8*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
                % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                % in the range, we may consider 20% or 30% reduction. Scatter correction saves this issue.
                Coordinates(counter,:)=[x_Tx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
                    y_Tx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                    z_Tx + a_c(counter)*sind(theta_av(counter))] ;
            end
            
        end
        
    elseif env==2 % Outdoors
        % Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
        Coordinates=zeros(C,3);          % for Clusters
        Coordinates2=zeros(sum(S),3);    % for Scatterers
        for counter=1:C
            loop=1;
            Coordinates(counter,:)=[x_Tx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
                y_Tx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                z_Tx + a_c(counter)*sind(theta_av(counter))] ;
            while  Coordinates(counter,3)<0    % only underground clusters will be ignored
                a_c(counter)=    0.8*a_c(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
                % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                % in the range, we may consider 20% or 30% reduction.
                Coordinates(counter,:)=[x_Tx + a_c(counter)*cosd(theta_av(counter))*cosd(phi_av(counter)),...
                    y_Tx - a_c(counter)*cosd(theta_av(counter))*sind(phi_av(counter)),...
                    z_Tx + a_c(counter)*sind(theta_av(counter))] ;
            end
        end
    end
    a_c_rep=[];
    for counter3=1:C
        a_c_rep=[a_c_rep,repmat(a_c(counter3),1,S(counter3))];
    end
    for counter2=1:sum(S)
        Coordinates2(counter2,:)=[x_Tx + a_c_rep(counter2)*cosd(theta_Tx(counter2))*cosd(phi_Tx(counter2)),...
            y_Tx - a_c_rep(counter2)*cosd(theta_Tx(counter2))*sind(phi_Tx(counter2)),...
            z_Tx + a_c_rep(counter2)*sind(theta_Tx(counter2))] ;
        
    end
    
    % Correction on Scatters
    % You may ignore the scatterers outside the walls for Enviroment 1 (Indoors) or underground for Enviroment 2 (Outdoors)
    
    if env==1
        ignore=[];
        
        for counter2=1:sum(S)
            if Coordinates2(counter2,3)>dim(3) || Coordinates2(counter2,3)<0 ||  Coordinates2(counter2,2)>dim(2) ||  Coordinates2(counter2,2)<0  ||  Coordinates2(counter2,1)>dim(1) ||  Coordinates2(counter2,1)<0
                ignore=[ignore,counter2];   % contains the indices of scatterers that can be ignored
            end
        end
        
        % updated indices
        indices=setdiff(1:sum(S),ignore);    % the set of active scatterer indices
        M_new=length(indices);               % number of IOs inside the room or above ground
        
        
        
    elseif env==2
        
        ignore=[];
        
        for counter2=1:sum(S)
            if  Coordinates2(counter2,3)<0   % only underground scatterers
                ignore=[ignore,counter2];   % contains the indices of scatterers that can be ignored
            end
        end
        
        % updated indices
        indices=setdiff(1:sum(S),ignore);    % the set of active scatterer indices
        M_new=length(indices);               % number of IOs inside the room or above ground
        
    end
    
    % if you want to revert back Version 1.1 from Version 1.2, simply set
    % indices=1:sum(S);
    % M_new=sum(S);
    
    % Neccessary Loop to have at least one scatterer
    if M_new>0 % if M_new==0 --> all scatters are outside
        break  % break generate=1:100 % if M_new >0 we are OK (at least one scatter)
    end
    
end  % for generate=1:100



% Generated h_LOS can be used for both envs.


phi_RIS=zeros(1,M_new);
theta_RIS=zeros(1,M_new);
b_cs=zeros(1,M_new);
d_cs=zeros(1,M_new);


normal = 0;
for counter2=indices
    normal = normal + 1;
    b_cs(counter2)=norm(RIS_xyz-Coordinates2(counter2,:));   % Distance between Scatterer and RIS
    d_cs(normal)=a_c_rep(counter2)+b_cs(counter2);         % Total distance Tx-Scatterer-RIS
    
    I_phi=sign(x_RIS-Coordinates2(counter2,1));
    phi_RIS(normal)  = I_phi* atand ( abs( x_RIS-Coordinates2(counter2,1)) / abs(y_RIS-Coordinates2(counter2,2)) );
    
    I_theta=sign(Coordinates2(counter2,3)-z_RIS);
    theta_RIS(normal)=I_theta * asind ( abs (z_RIS-Coordinates2(counter2,3) )/ b_cs(counter2) );
    
end
if I_LOS == 1
    phi_RIS(1) = phi_RIS_LOS;
    theta_RIS(1) = theta_RIS_LOS;
end



h_ncs = zeros(1,2,N,M_new);
Pc = 1 ;
GeT = [1;1]./sqrt(2);
beta=zeros(1,M_new); % to be reused for shared clusters (env 1)
shadow=beta;          % to be reused for shared clusters (env 1) - Version 1.4
alfa_cs_siso = zeros(2,2,M_new);


for counter2=1:M_new %toplam cs uzerinden los dahil
    counter3=1; %toplam n uzerinden donuyor
    
    
    if I_LOS == 1 && counter2 ==1
        gama = 1;
        I = 1;
        Lcs = L_LOS;
        alfa_cs = [exp(1i*rand*2*pi) 0; 0 exp(1i*rand*2*pi)];
        beta(counter2) = 1;
        shadow(counter2) = hsiso_val;
    else
        gama = sqrt(1/M_new);
        beta(counter2) = ((randn+1i*randn)./sqrt(2));
        I = 1;
        X_sigma=randn*sigma_NLOS;
        shadow(counter2)=X_sigma;
        Lcs_dB=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs(counter2))- X_sigma;
        Lcs=10^(Lcs_dB/10);
        X = db2pow(xpr_std_NLOS)*randn + db2pow((xpr_mu_NLOS));
        Kappa = 10^(X/10);
        alfa_cs = [exp(1i*rand*2*pi) sqrt(1/Kappa)*exp(1i*rand*2*pi); sqrt(1/Kappa)*exp(1i*rand*2*pi) exp(1i*rand*2*pi)];
        alfa_cs_siso(:,:,counter2)=alfa_cs; %for h siso
    end
    
    for x=0:sqrt(N)-1
        for y=0:sqrt(N)-1
            h_ncs(:,:,counter3,counter2)= gama * I * sqrt(Lcs) * abs(beta(counter2)) * transpose(GeT) * alfa_cs * array_resp(theta_RIS(counter2),phi_RIS(counter2),x,y);
            counter3=counter3+1;
        end
    end
end

%% g kanali
if env==1   % PURE LOS CHANNEL
    
    d_RIS_R=norm(RIS_xyz-Rx_xyz);
    
    I_theta=sign(z_Rx - z_RIS);
    theta_RIS_LOS_2=I_theta * asind( abs(z_Rx-z_RIS)/d_RIS_R );
    I_phi=sign(x_RIS - x_Rx);
    phi_RIS_LOS_2=I_phi * atand( abs(x_Rx-x_RIS)/ abs(y_Rx-y_RIS) );
    
    array_LOS_2=zeros(1,N);
    
    counter3=1;
    for x=0:sqrt(N)-1
        for y=0:sqrt(N)-1
            array_LOS_2(counter3)=exp(1i*k*dis*(x*sind(theta_RIS_LOS_2) + y*sind(phi_RIS_LOS_2)*cosd(theta_RIS_LOS_2) )) ;
            counter3=counter3+1;
        end
    end
    
    L_dB_LOS_2=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_RIS_R)- randn*sigma_LOS;
    L_LOS_2=10^(L_dB_LOS_2/10);
    
    alfa_g_LOS = [exp(1i*rand*2*pi) 0 ; 0 exp(1i*rand*2*pi)];
    GeR = [1;1]./sqrt(2); %Omnidirectional Pattern assumed. Other antenna patterns can be used by considering Tx-RIS angle.
    g_LOS = zeros(2,1,N);
    
    for n=1:N
        g_LOS(:,:,n)=sqrt(L_LOS_2)*alfa_g_LOS*GeR*array_LOS_2(n);
    end
    
    g_n = g_LOS;
    
elseif env==2
    
    d_RIS_R=norm(RIS_xyz-Rx_xyz);
    
    p_LOS_2=min([20/d_RIS_R,1])*(1-exp(-d_RIS_R/39)) + exp(-d_RIS_R/39);
    
    I_LOS_2=randsrc(1,1,[1,0;p_LOS_2 1-p_LOS_2]);
    
    if I_LOS_2==1
        
        I_phi=sign(x_RIS-x_Rx);
        phi_RIS_LOS_2 = I_phi* atand ( abs( x_RIS-x_Rx) / abs(y_RIS-y_Rx) );
        I_theta=sign(z_Rx-z_RIS);
        theta_RIS_LOS_2=I_theta * asind ( abs (z_RIS-z_Rx )/ d_RIS_R );
        
        L_dB_LOS_2=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_RIS_R)- randn*sigma_LOS;
        L_LOS_2=10^(L_dB_LOS_2/10);
    else
        g_LOS=0;
        phi_RIS_LOS_2 = 0;
        theta_RIS_LOS_2 = 0;
    end
    
    for generate2=1:100  % To ensure that at least one scatterer exist
        
        if I_LOS_2 == 1
            C_2=max([1,poissrnd(lambda_p)])+1;  % Poisson distributed
            S_2=randi(30,1,C_2);
            S_2(1) = 1;
        else
            C_2=max([1,poissrnd(lambda_p)]);  % Poisson distributed
            S_2=randi(30,1,C_2);
        end
        
        phi_Tx_2=[ ];
        theta_Tx_2=[ ];
        phi_av_2=zeros(1,C_2);
        theta_av_2=zeros(1,C_2);
        
        for counter=1:C_2
            
            phi_av_2(counter)  = rand*90-45;   % mean azimuth (reduced to ensure that all scatters are within the field of view)
            theta_av_2(counter)= rand*90-45; % mean elevation  (needs update & most of the clusters are above ceiling)
            % for outdoors this might not be a big concern
            
            % cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
            phi_Tx_2   = [phi_Tx_2,  log(rand(1,S_2(counter))./rand(1,S_2(counter)))*sqrt(25/2) + phi_av_2(counter)];
            theta_Tx_2 = [theta_Tx_2,log(rand(1,S_2(counter))./rand(1,S_2(counter)))*sqrt(25/2) + theta_av_2(counter)];
            
        end
        
        
        
        % Cluster Distances
        a_c_2=1+rand(1,C_2)*(d_RIS_R-1);    % Cluster distances uniform [1,d_RIS_R] % can be modified later
        
        
        % Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
        
        Coordinates_2=zeros(C_2,3);          % for Clusters
        Coordinates2_2=zeros(sum(S_2),3);    % for Scatterers
        
        for counter=1:C_2
            loop=1;
            
            Coordinates_2(counter,:)=[x_RIS - a_c_2(counter)*cosd(theta_av_2(counter))*sind(phi_av_2(counter)),...
                y_RIS - a_c_2(counter)*cosd(theta_av_2(counter))*cosd(phi_av_2(counter)),...
                z_RIS + a_c_2(counter)*sind(theta_av_2(counter))];
            
            while  Coordinates_2(counter,3)<0    % only underground clusters will be ignored
                a_c_2(counter)=    0.8*a_c_2(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
                Coordinates_2(counter,:)=[x_RIS - a_c_2(counter)*cosd(theta_av_2(counter))*sind(phi_av_2(counter)),...
                    y_RIS - a_c_2(counter)*cosd(theta_av_2(counter))*cosd(phi_av_2(counter)),...
                    z_RIS + a_c_2(counter)*sind(theta_av_2(counter))];
                
            end
            
        end
        
        
        a_c_rep_2=[];
        for counter3=1:C_2
            a_c_rep_2=[a_c_rep_2,repmat(a_c_2(counter3),1,S_2(counter3))];
        end
        
        
        for counter2=1:sum(S_2)
            Coordinates2_2(counter2,:)=[x_RIS - a_c_rep_2(counter2)*cosd(theta_Tx_2(counter2))*sind(phi_Tx_2(counter2)),...
                y_RIS - a_c_rep_2(counter2)*cosd(theta_Tx_2(counter2))*cosd(phi_Tx_2(counter2)),...
                z_RIS + a_c_rep_2(counter2)*sind(theta_Tx_2(counter2))];
        end
        
        ignore_2=[];
        
        for counter2=1:sum(S_2)
            if  Coordinates2_2(counter2,3)<0   % only underground scatterers
                ignore_2=[ignore_2,counter2];   % contains the indices of scatterers that can be ignored
            end
        end
        
        % updated indices
        indices_2=setdiff(1:sum(S_2),ignore_2);    % the set of active scatterer indices
        M_new_2=length(indices_2);               % number of IOs inside the room or above ground
        
        
        % Neccessary Loop to have at least one scatterer
        if M_new_2>0 % all scatters are outside
            break  % break generate=1:100
        end
        
    end  % for generate=1:100
    
    
    % Calculate Link Lengths (Step 3) and Path Loss (Step 5)
    b_cs_2=zeros(1,M_new_2);
    d_cs_2=zeros(1,M_new_2);
    
    
    normal = 0;
    for counter2=indices_2
        normal = normal + 1;
        b_cs_2(counter2)=norm(Rx_xyz-Coordinates2_2(counter2,:));   % Distance between Scatterer and Rx
        d_cs_2(normal)=a_c_rep_2(counter2)+b_cs_2(counter2);         % Total distance RIS-Scatterer-Rx
        
    end
    
    if I_LOS_2 == 1
        phi_Tx_2(1) = phi_RIS_LOS_2;
        theta_Tx_2(1) = theta_RIS_LOS_2;
        
    end
    
    g_nkl = zeros(2,1,N,M_new_2);
    GeR = [1;1]./sqrt(2);
    
    
    for counter2=1:M_new_2
        counter3=1;
        
        if counter2 == 1 && I_LOS_2 == 1 %LOS CONDITION
            
            Ig = 1;
            Lcs_2 = L_LOS_2;
            alfa_cs_2 = [exp(1i*rand*2*pi) 0; 0 exp(1i*rand*2*pi)];
            gama2 = 1;
            beta_2 =1;
            
            
        else %NLOS CONDITION
            Ig = 1;
            gama2 = sqrt(1/M_new_2);
            beta_2=((randn+1i*randn)./sqrt(2));
            X_sigma_2=randn*sigma_NLOS;
            Lcs_dB_2=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_2(counter2))- X_sigma_2;
            Lcs_2=10^(Lcs_dB_2/10);
            X_2 = db2pow(xpr_std_NLOS)*randn + db2pow((xpr_mu_NLOS));
            Kappa_2 = 10^(X_2/10);
            alfa_cs_2 = [exp(1i*rand*2*pi) sqrt(1/Kappa_2)*exp(1i*rand*2*pi); sqrt(1/Kappa_2)*exp(1i*rand*2*pi) exp(1i*rand*2*pi)];
        end
        
        for x=0:sqrt(N)-1
            for y=0:sqrt(N)-1
                g_nkl(:,:,counter3,counter2)=gama2*Ig*sqrt(Lcs_2)*abs(beta_2)*alfa_cs_2*GeR*array_resp(theta_Tx_2(counter2),phi_Tx_2(counter2),x,y);
                counter3=counter3+1;
                
            end
        end
        
    end
    
end

%% STEP 8
% Generation of h_SISO
if siso == 1
    if env==1  % WITH SHARED CLUSTERS AND COMMON LOS
        
        d_T_R=norm(Tx_xyz-Rx_xyz);
        
        d_cs_tilde=zeros(1,M_new);
        h_SISO_NLOS=0;
        
        
        
        if I_LOS == 1
            indices(1) = [];
        end
        
        sayac = 1;
        for counter2=indices
            
            % due to shared clusters d_cs_tilde ~ d_cs
            d_cs_tilde(counter2) = a_c_rep(counter2) + norm(Coordinates2(counter2,:)- Rx_xyz);
            
            %    X_sigma=randn*sigma; SHARED CLUSTERS (same Shadow Factor - Version 1.4 Updated)
            
            Lcs_dB_SISO=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_tilde(counter2))- shadow(sayac);
            
            Lcs_SISO=10^(Lcs_dB_SISO/10);
            
            % We consider the same complex path gain (small scale fading) with an excess phase (similar to array response)
            eta=k* ( norm(Coordinates2(counter2,:)- RIS_xyz) -  norm(Coordinates2(counter2,:)- Rx_xyz));
            
            h_SISO_NLOS = h_SISO_NLOS + abs(beta(sayac))*exp(1i*eta)*sqrt(Lcs_SISO)*transpose(GeT)*alfa_cs_siso(:,:,sayac)*GeR;
            sayac = sayac + 1;
        end
        
        h_SISO_NLOS=h_SISO_NLOS.*sqrt(1/M_new);  % normalization
        
        if z_RIS >= z_Tx
            
            % % Include LOS component (Version 1.4)
            %     % InH LOS Probability
            if d_T_R<= 1.2
                p_LOS_3=1;
            elseif 1.2<d_T_R && d_T_R<6.5
                p_LOS_3=exp(-(d_T_R-1.2)/4.7);
            else
                p_LOS_3=0.32*exp(-(d_T_R-6.5)/32.6);
            end
            
            I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);
            
            % Do not recalculate, if T-RIS has LOS, we might have LOS for h_SISO as well (for ground level RIS)
            % If we have LOS for Tx-RIS, we have LOS for Tx-Rx
        elseif z_RIS < z_Tx  % RIS in the ground level
            
            I_LOS_3=I_LOS;
            
        end
        
        if I_LOS_3==1
            alfa_cs_3 = [exp(1i*rand*2*pi) 0; 0 exp(1i*rand*2*pi)];
            L_SISO_LOS_dB=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- randn*sigma_LOS;
            L_SISO_LOS=10^(L_SISO_LOS_dB/10);
            h_SISO_LOS= sqrt(L_SISO_LOS)*transpose(GeT)*alfa_cs_3*GeR;
        else
            h_SISO_LOS=0;
        end
        
        h_SISO=h_SISO_NLOS + h_SISO_LOS; % include LOS Component if I_LOS_3==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif env==2 % WITH NEW CLUSTERS
        
        % Generate New Clusters/Scatterers (Outdoors)
        d_T_R=norm(Tx_xyz-Rx_xyz);
        % Calculate LOS Component (Version 1.4)
        % UMi LOS Probability
        p_LOS_3=min([20/d_T_R,1])*(1-exp(-d_T_R/39)) + exp(-d_T_R/39);
        I_LOS_3=randsrc(1,1,[1,0;p_LOS_3 1-p_LOS_3]);
        if I_LOS_3==1
            alfa_cs_3 = [exp(1i*rand*2*pi) 0; 0 exp(1i*rand*2*pi)];
            L_SISO_LOS_dB=-20*log10(4*pi/lambda) - 10*n_LOS*(1+b_LOS*((Frequency-f0)/f0))*log10(d_T_R)- randn*sigma_LOS;
            L_SISO_LOS=10^(L_SISO_LOS_dB/10);
            h_SISO_LOS= sqrt(L_SISO_LOS)*transpose(GeT)*alfa_cs_3*GeR;   % LOS Component for the Tx-Rx link
        else
            h_SISO_LOS=0;
        end
        % Generate New Clusters/Scatters (Step 2) - Ensure that All Scattters are above ground
        for generate3=1:100  % To ensure that at least one scatterer exist
            % lambda_p was defined before for 28/73 GHz
            C_3=max([1,poissrnd(lambda_p)]);  % Poisson distributed
            % Number of Sub-rays per Cluster
            S_3=randi(30,1,C_3);
            % Azimuth/Elevation Departure Angles
            phi_Tx_3=[ ];
            theta_Tx_3=[ ];
            phi_av_3=zeros(1,C_3);
            theta_av_3=zeros(1,C_3);
            for counter=1:C_3
                phi_av_3(counter)  = rand*180-90;   % mean azimuth
                theta_av_3(counter)= rand*90-45; % mean elevation  (needs update & most of the clusters are above ceiling)
                % for outdoors this might not be a big concern
                % cluster angles: First S(1) belongs to Cluster 1, Next S(2) belongs to Cluster 2....
                phi_Tx_3   = [phi_Tx_3,  log(rand(1,S_3(counter))./rand(1,S_3(counter)))*sqrt(25/2) + phi_av_3(counter)];
                theta_Tx_3 = [theta_Tx_3,log(rand(1,S_3(counter))./rand(1,S_3(counter)))*sqrt(25/2) + theta_av_3(counter)];
            end
            % Cluster Distances
            a_c_3=1+rand(1,C_3)*(d_T_R-1);    % Cluster distances uniform [1,d_T_R]
            % Reducing Distances for Outside Clusters - Check Cluster Coordinates with Mean Angles
            Coordinates_3=zeros(C_3,3);          % for Clusters
            Coordinates2_3=zeros(sum(S_3),3);    % for Scatterers
            for counter=1:C_3
                loop=1;
                Coordinates_3(counter,:)=[x_Tx + a_c_3(counter)*cosd(theta_av_3(counter))*cosd(phi_av_3(counter)),...
                    y_Tx - a_c_3(counter)*cosd(theta_av_3(counter))*sind(phi_av_3(counter)),...
                    z_Tx + a_c_3(counter)*sind(theta_av_3(counter))] ;
                while  Coordinates_3(counter,3)<0    % only underground clusters will be ignored
                    a_c_3(counter)=    0.8*a_c_3(counter)  ;     % reduce distance 10%-20% if coordinates are not met!
                    % Note: While the above 10% reduction ensures that all clusters are in the range, to ensure that most of the scatterers are
                    % in the range, we may consider 20% or 30% reduction.
                    
                    Coordinates_3(counter,:)=[x_Tx + a_c_3(counter)*cosd(theta_av_3(counter))*cosd(phi_av_3(counter)),...
                        y_Tx - a_c_3(counter)*cosd(theta_av_3(counter))*sind(phi_av_3(counter)),...
                        z_Tx + a_c_3(counter)*sind(theta_av_3(counter))] ;
                    %    loop=loop+1
                end
                %     % Plot Clusters (comment figures in mass generation of channels (with repeat))
                % plot3(Coordinates_3(counter,1),Coordinates_3(counter,2),Coordinates_3(counter,3),'rdiamond')
                % hold on;
            end
            % Plot Scatterers
            a_c_rep_3=[];
            for counter3=1:C_3
                a_c_rep_3=[a_c_rep_3,repmat(a_c_3(counter3),1,S_3(counter3))];
            end
            for counter2=1:sum(S_3)
                Coordinates2_3(counter2,:)=[x_Tx + a_c_rep_3(counter2)*cosd(theta_Tx_3(counter2))*cosd(phi_Tx_3(counter2)),...
                    y_Tx - a_c_rep_3(counter2)*cosd(theta_Tx_3(counter2))*sind(phi_Tx_3(counter2)),...
                    z_Tx + a_c_rep_3(counter2)*sind(theta_Tx_3(counter2))] ;
                %       % comment figures in mass generation of channels (with repeat)
                %     plot3(Coordinates2_3(counter2,1),Coordinates2_3(counter2,2),Coordinates2_3(counter2,3),'r.')
                %   hold on;
            end
            % Correction on Scatters
            ignore_3=[];
            for counter2=1:sum(S_3)
                if  Coordinates2_3(counter2,3)<0   % only underground scatterers
                    ignore_3=[ignore_3,counter2];   % contains the indices of scatterers that can be ignored
                end
            end
            % updated indices
            indices_3=setdiff(1:sum(S_3),ignore_3);    % the set of active scatterer indices
            M_new_3=length(indices_3);               % number of IOs inside the room or above ground
            % Neccessary Loop to have at least one scatterer
            if M_new_3>0 % all scatters are outside
                break  % break generate=1:100
            end
        end  % for generate=1:100
        % Calculate Link Lengths and Path Loss
        b_cs_3=zeros(1,sum(S_3));
        d_cs_3=zeros(1,sum(S_3));
        for counter2=indices_3
            b_cs_3(counter2)=norm(Tx_xyz-Coordinates2_3(counter2,:));      % Distance between Scatterer and Tx
            d_cs_3(counter2)=a_c_rep_3(counter2)+b_cs_3(counter2);         % Total distance Tx-Scatterer-Rx
        end
        % Generate h_SISO
        h_SISO_NLOS=0;
        for counter2=indices_3
            X_sigma_3=randn*sigma_NLOS;
            Lcs_dB_SISO_3=-20*log10(4*pi/lambda) - 10*n_NLOS*(1+b_NLOS*((Frequency-f0)/f0))*log10(d_cs_3(counter2))- X_sigma_3;
            Lcs_SISO_3=10^(Lcs_dB_SISO_3/10);
            X_3 = db2pow(xpr_std_NLOS)*randn + db2pow((xpr_mu_NLOS));
            Kappa_3 = 10^(X_3/10);
            alfa_cs_3 = [exp(1i*rand*2*pi) sqrt(1/Kappa_3)*exp(1i*rand*2*pi); sqrt(1/Kappa_3)*exp(1i*rand*2*pi) exp(1i*rand*2*pi)];
            h_SISO_NLOS = h_SISO_NLOS + ((randn+1i*randn)/sqrt(2))*sqrt(Lcs_SISO_3)*transpose(GeT)*alfa_cs_3*GeR;
        end
        h_SISO_NLOS = h_SISO_NLOS.*sqrt(1/M_new_3);  % normalization
        h_SISO = h_SISO_NLOS + h_SISO_LOS;    % h_SISO_LOS=0  if there is no LOS
    end
    
else
    h_SISO = 0;
end

if env == 2
    H = 0;
    for nn = 1:N
        for cs=1:M_new
            for kl = 1:M_new_2
                %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_nkl(:,:,nn,kl))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn,kl)));
                H_slack = abs(h_ncs(:,:,nn,cs))* abs(g_nkl(:,:,nn,kl))*exp(1i*angle(h_SISO)) ;
                H = H + H_slack;
            end
        end
    end
    H = H + h_SISO;
else
    H = 0;
    for nn = 1:N
        for cs=1:M_new
            %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_n(:,:,nn))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn)));
            H_slack = abs(h_ncs(:,:,nn,cs))* abs(g_n(:,:,nn))*exp(1i*angle(h_SISO)) ;
            H = H + H_slack;
        end
    end
    H = H + h_SISO;
end
end