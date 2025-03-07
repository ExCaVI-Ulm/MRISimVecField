function [M_propG] = PropMatrixBconstSum_block(settings,indices_start, indices_end, B_t_woRF, time_tot, M_prev, j)
%PropMatrixBconstSum_block blockwise propagation with matrix-based
%approach WITHOUT relaxation contributions
%   Input:  -settings struct
%           -indices_start: starting indices of the individual blocks
%           -indices_end: indices marking the ends of the individual blocks
%           -B_t_woRF:  Time course of Magnetic fields without excitation contributions
%           -time_tot: time vector 
%           -M_prev:  magnetization prior to this block of propagation
%           -j: current number of block being processed 
%   Output:
%           -M_propG: propagated magnetization vector after application of
%           constant magnetic vector field

    matrixsize_signal = settings.signal.matrixsize_signal;

    M_propG = zeros(indices_end(j) - indices_start(j)+1, matrixsize_signal^2,3);

    gamma = settings.general.gamma;
    sampling_freq = settings.signal.sampling_freq;
    parfor ll = 1:matrixsize_signal^2

        B_t_prop = squeeze(B_t_woRF(1,ll,:));
        u_B_t_prop = B_t_prop / norm(B_t_prop);

        x1_r = 1/sqrt(2*(u_B_t_prop(1)^2 + u_B_t_prop(3)^2))*[u_B_t_prop(1)*u_B_t_prop(2) - 1i*u_B_t_prop(3) ; - (u_B_t_prop(1)^2 + u_B_t_prop(3)^2) ; u_B_t_prop(2)*u_B_t_prop(3) + 1i*u_B_t_prop(1)];
        x2_r = u_B_t_prop;
        x3_r = conj(x1_r);

        B_tp_mag_r_t = norm(B_t_prop);
        lambda3_r_t = 1i*gamma*(B_tp_mag_r_t)*(1:(indices_end(j)-indices_start(j)+1)).'*1/sampling_freq;

        lambda2 = 0;
        lambda1_r_t = conj(lambda3_r_t);

        M_t0 = M_prev;
        v_1 = exp(lambda1_r_t)*x1_r'*squeeze(M_t0(ll,:)).'*x1_r(1) + exp(lambda2)*x2_r'*squeeze(M_t0(ll,:)).'*x2_r(1) + exp(lambda3_r_t)*x3_r'*squeeze(M_t0(ll,:)).'*x3_r(1);
        v_2 = exp(lambda1_r_t)*x1_r'*squeeze(M_t0(ll,:)).'*x1_r(2) + exp(lambda2)*x2_r'*squeeze(M_t0(ll,:)).'*x2_r(2) + exp(lambda3_r_t)*x3_r'*squeeze(M_t0(ll,:)).'*x3_r(2);
        v_3 = exp(lambda1_r_t)*x1_r'*squeeze(M_t0(ll,:)).'*x1_r(3) + exp(lambda2)*x2_r'*squeeze(M_t0(ll,:)).'*x2_r(3) + exp(lambda3_r_t)*x3_r'*squeeze(M_t0(ll,:)).'*x3_r(3);
        v = [v_1,v_2, v_3];
        M_propG(:,ll,:) = v; 
    end
end

