function [M_prop] = ExcMatrixCalc_Reco(settings, B_t_woRF, indices_start, indices_end, B1TX_spatial, blocks, j, Minit, reco_flag)
%ExcMatrixCalc_Reco Matrix based propagation for an RF block for reconstruction
%   Input:  -settings struct
%           -B_t_woRF:  Time course of Magnetic fields without excitation contributions
%           -indices_start: starting indices of the individual blocks
%           -indices_end: indices marking the ends of the individual blocks
%           -B1TX_spatial: spatially varying contributions of TX B1
%           -blocks: cell with all blocks of the simulated sequence
%           -j: current block
%           -Minit: initial magnetization
%           -reco_flag: boolean which enables to switch between function
%           usage for reconstruction or simulation
%   Output:
%           -M_prop: propagated Magnetization vector
    if reco_flag
        matrixsize_signal = settings.reco.matrixsize_reco^2;
    else
        matrixsize_signal = settings.signal.matrixsize_signal^2;
    end

    B1TX_spatial = reshape(B1TX_spatial,3,[]);

    u_B = zeros(matrixsize_signal,3);
    rot_axis = zeros(matrixsize_signal,3);
    u_B1_trafo = zeros(matrixsize_signal,3);
    M_exc = zeros(matrixsize_signal,3,settings.trajectory.N_PhaseEnc);
    w = settings.TX.w;
    for uuu = 1:settings.trajectory.N_PhaseEnc
    for ll = 1:matrixsize_signal
        %transformation into tilted B0-frame
        u_B(ll,:) = B_t_woRF(ll,:,uuu) / norm(squeeze(B_t_woRF(ll,:,uuu)));    %calculate unit vector

        %rot axis
        rot_axis(ll,:) = cross(squeeze(B_t_woRF(ll,:,uuu)), [0;0;1]);
        if norm(rot_axis(ll,:)) ~= 0
            rot_axis(ll,:) = rot_axis(ll,:) / norm(cross(squeeze(B_t_woRF(ll,:,uuu)), [0;0;1]));
        end
        theta = acos(squeeze(B_t_woRF(ll,3,uuu)) / norm(squeeze(B_t_woRF(ll,:,uuu))));
        F_mat = [ 0, -rot_axis(ll,3), rot_axis(ll,2); ...
                rot_axis(ll,3), 0, -rot_axis(ll,1); ...
                -rot_axis(ll,2), rot_axis(ll,1), 0];
        R_CoSy_theta = eye(3) + sin(theta)* F_mat + (1-cos(theta))*F_mat*F_mat;
        R_CoSy_Mtheta = eye(3) - sin(theta)* F_mat + (1-cos(theta))*F_mat*F_mat;
    
        u_B1_trafo(ll,:) = R_CoSy_theta*(squeeze(B1TX_spatial(:,ll)) / norm(squeeze(B1TX_spatial(:,ll))));%unit vector B1
        u_B1_trafo(ll,3) = 0;%
        %excitation
        a_rot_exc = [cos(blocks{j}.initPhase);sin(blocks{j}.initPhase);0] + (norm(squeeze(B_t_woRF(ll,:,uuu)))*settings.general.gamma - w(1))/settings.general.gamma * [0;0;1];
        
        alpha = norm(u_B1_trafo(ll,:)) *blocks{j}.flipangle;
    
        A_mat = [0, -a_rot_exc(3), a_rot_exc(2) ; a_rot_exc(3), 0, -a_rot_exc(1); -a_rot_exc(2) , a_rot_exc(1), 0];
        R_exc_alpha = eye(3) - sin(alpha)*A_mat + (1-cos(alpha))*A_mat*A_mat;
        
        %rotating frame
        phase =  -abs(w(1))*blocks{j}.dur;
        R_rotframe = [cos(phase), -sin(phase), 0 ; sin(phase), cos(phase) 0; 0, 0, 1];
        R_rotframeB = [cos(-phase), -sin(-phase), 0 ; sin(-phase), cos(-phase) 0; 0, 0, 1];
        
        M_exc(ll,:,uuu) = R_CoSy_Mtheta*R_rotframe*R_exc_alpha*R_rotframeB*R_CoSy_theta*Minit(ll,:,uuu).';
    end
    end

    M_prop = permute(repmat(M_exc, [1 1 1 (indices_end(j)-indices_start(j)+1)]), [4 1 2 3]);
end