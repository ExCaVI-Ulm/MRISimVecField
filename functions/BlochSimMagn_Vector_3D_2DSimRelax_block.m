function [M] = BlochSimMagn_Vector_3D_2DSimRelax_block(settings, time_tot, B_t, sampleS, Minit)
%BlochSimMagn_Vector_3D_2DSimRelax_block Bloch Simulation with relaxation and with arbitrary magnetic vector field B_t and initial magnetization Minit
%   Input:  -settings struct
%           -time vector time_tot
%           -B_t: Magnetic vector field (possibly time dependent<-> no reconstruction for that) during this block
%           -sampleS: struct with informations about the sample - M0, T1,T2
%           -Minit: initial magnetization
%   Output:
%           -M: Magnetization vector after simulation

    gamma_loc = settings.general.gamma;
    matrixsize = settings.signal.matrixsize_signal;
    M = zeros(length(time_tot), settings.signal.matrixsize_signal^2, 3);
        
    %dt_Bloch = time_tot(end)/length(time_tot);
    dt_Bloch = time_tot(2)-time_tot(1);

    T1 = reshape(sampleS.T1, [],1);
    T2 = reshape(sampleS.T2, [],1);
    
    Minit_straight = reshape(Minit, [], 3);
    NF = length(time_tot);
    for j =1:(matrixsize)^2%par
        % if (Minit(j,1) == 0 && Minit(j,2) == 0 && Minit(j,3) == 0)
        %     continue;
        % end
        M_sim = zeros(3,2);
        M_sim(:,1) = squeeze(Minit_straight(j,:)).';    

        for u=1:NF
            Btot_z = B_t(u,j,3);
            Btot_y = B_t(u,j,2);
            Btot_x = B_t(u,j,1);
             A = [0,  gamma_loc*Btot_z,  -gamma_loc*Btot_y; ...        
                  -gamma_loc*Btot_z,   0,  gamma_loc*Btot_x; ...
                   gamma_loc*Btot_y,   -gamma_loc*Btot_x,     0];
                    
            P_B = [Btot_x; Btot_y; Btot_z]*[Btot_x, Btot_y, Btot_z] / (dot([Btot_x; Btot_y; Btot_z],[Btot_x; Btot_y; Btot_z]));
            M_sim(:,2) = M_sim(:,1) + ((A- eye(3)/T2(j) + (1/T2(j) - 1/T1(j))*P_B )*M_sim(:,1))*dt_Bloch + squeeze(Minit_straight(j,:)).' / (T1(j));
            M_sim(:,1) = M_sim(:,2); 
            M(u,j,:) = M_sim(1:3,end);
        end
        
    end

end

