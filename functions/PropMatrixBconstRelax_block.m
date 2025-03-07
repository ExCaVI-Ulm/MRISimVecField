function [M_propG] = PropMatrixBconstRelax_block(settings,indices_start, indices_end, B_t_woRF, time_tot, M_prev, sampleS, Minit, j)
%PropMatrixBconstRelax_block blockwise propagation with matrix-based
%approach using relaxation contributions
%   Input:  -settings struct
%           -indices_start: starting indices of the individual blocks
%           -indices_end: indices marking the ends of the individual blocks
%           -B_t_woRF:  Time course of Magnetic fields without excitation contributions
%           -time_tot: time vector 
%           -M_prev:  magnetization prior to this block of propagation
%           -sampleS: sample struct with information about magnetization,
%           T1, T2,...
%           -Minit: initial magnetization vector
%           -j: current number of block being processed 
%   Output:
%           -M_propG: propagated magnetization vector after application of
%           constant magnetic vector field

    matrixsize_signal = settings.signal.matrixsize_signal;

    T2 = reshape(sampleS.T2,[],1);
    T1 = reshape(sampleS.T1,[],1);

    Minit = reshape(Minit,[],3);

    M_propG = zeros(indices_end(j) - indices_start(j)+1, matrixsize_signal^2,3);

    for ll = 1:matrixsize_signal^2
        B_t_prop = squeeze(B_t_woRF(1,ll,:));
        A_prop = settings.general.gamma*...
                 [ 0, B_t_prop(3), -B_t_prop(2); ...
                  -B_t_prop(3), 0,  B_t_prop(1); ...
                   B_t_prop(2), -B_t_prop(1), 0];
    
        Projmatrix = B_t_prop * B_t_prop.' / dot(B_t_prop,B_t_prop);
        A_prop = A_prop -  1/T2(ll)* eye(size(A_prop)) + (1/T2(ll) - 1/T1(ll)) * Projmatrix;
        
        A_bar3 = [A_prop squeeze(Minit(ll,:)).'/T1(ll)];
        A_bar4 = [A_bar3; zeros(1,4)];

        t_prop = time_tot(indices_start(j):indices_end(j));
        if j ~= 1
            t_prop = t_prop - time_tot(indices_end(j-1));
        end
    
        M_t0 = M_prev;
        % for uu = 1:length(t_prop)
        %     M_propG(uu, ll,:) =expm(A_prop*t_prop(uu))*squeeze(M_t0(ll,:)).';
        % end
        %tic
        M_add_relax_contr = expmv(A_bar4, [0;0;0;1], t_prop).';
        M_propG(:,ll,:) = expmv(A_prop, squeeze(M_t0(ll,:)).', t_prop).' + M_add_relax_contr(1:3);
        %toc
        clearvars A_bar
    end

end