function H = channel_calculateVAR(C)

global N M_new M_new_2 theta_RIS status hncs gn h_SISO env

if status == "ideal"
    if env == 1
        H = 0;
        for nn = 1:N
            for cs=1:M_new
                %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_n(:,:,nn))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn)));
                H_slack = abs(hncs(:,:,nn,cs))* abs(gn(:,:,nn))*exp(1i*angle(h_SISO)) ;
                H = H + H_slack;
            end
        end
        H = H + h_SISO;
        H = (-abs(H)^2);
        
    elseif env == 2
        H = 0;
        for nn = 1:N
            for cs=1:M_new
                for kl = 1:M_new_2
                    %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_nkl(:,:,nn,kl))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn,kl)));
                    H_slack = abs(hncs(:,:,nn,cs))* abs(gn(:,:,nn,kl))*exp(1i*angle(h_SISO)) ;
                    H = H + H_slack;
                end
            end
        end
        H = H + h_SISO;
        H = (-abs(H)^2);
    else
        error("Please specify the environment");
    end
    
elseif status == "realistic"
    C= C*1e-10;
    if env == 1
        H = 0;
        for nn = 1:N
            for cs=1:M_new
                %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_n(:,:,nn))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn)));
                H_slack = hncs(:,:,nn,cs)*RIS_denemeVAR(theta_RIS(cs),C(nn))*gn(:,:,nn) ;
                H = H + H_slack;
            end
        end
        H = H + h_SISO;
        H = (-abs(H)^2);
        H = H*10^(16);
        
    elseif env == 2
        H = 0;
        for nn = 1:N
            for cs=1:M_new
                for kl = 1:M_new_2
                    %theta_ncs = transpose(conj(h_ncs(:,:,nn,cs)))*transpose(conj(g_nkl(:,:,nn,kl))) / (abs(h_ncs(:,:,nn,cs)) *abs(g_nkl(:,:,nn,kl)));
                    H_slack = hncs(:,:,nn,cs)*RIS_denemeVAR(theta_RIS(cs),C(nn))*gn(:,:,nn,kl) ;
                    H = H + H_slack;
                end
            end
        end
        H = H + h_SISO;
        H = (-abs(H)^2);
    else
        error("Please specify the environment");
    end
else
    error("Response can be ideal or realistic, please denote the case");
end

end