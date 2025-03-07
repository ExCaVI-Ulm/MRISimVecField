function [settings, B1vec_out, RF_blocks_ind] = TimeDependenceB12D(settings, blocks, indices_start, indices_end, B1TX_spatial, time_tot, sampling_freq)
%TimeDependenceB12D Calculation of time dependence of the B1 field for TX
%   Input:  -settings struct
%           -blocks: all blocks being simulated
%           -indices_start: starting indices of the individual blocks
%           -indices_end: indices marking the ends of the individual blocks
%           -B1TX_spatial: spatially varying contributions of TX B1
%           -time_tot: time vector
%           -sampling_freq: sampling frequency
%   Output:
%           -settings: settings struct
%           -B1vec_out: B1 as used for TX (now time dependent)
%           -RF_blocks_ind: indices of all RF blocks as defined in the
%           beginning of "SimMain"

    num_blocks = settings.general.num_blocks;
    matrixsize_signal = settings.signal.matrixsize_signal;
    matrixsize_reco = settings.reco.matrixsize_reco;

    for j = 1:num_blocks
        eval(['block',num2str(j),' = blocks{', num2str(j), '};'])
    end

    count = 0;
    for jj = 1:num_blocks
        if blocks{jj}.RF_block                      
            count = count+1;
            RF_blocks_ind(jj) =  jj;
        end
    end
    settings.trajectory.nr_RFblocks =count;
    RF_blocks_ind = nonzeros(RF_blocks_ind).';
    time_length = zeros(length(RF_blocks_ind),1);
    w = settings.general.gamma*settings.general.B0;
    settings.TX.w = w;
    
    for j=1:length(RF_blocks_ind)
        time_length(j) = blocks{RF_blocks_ind(j)}.dur;  
        
        eval(['timeRF_block', num2str(RF_blocks_ind(j)), '= time_tot(indices_start(RF_blocks_ind(', num2str(j), ')) : indices_end(RF_blocks_ind(', num2str(j), '))) - time_tot(indices_start(RF_blocks_ind(', num2str(j), '))) + 1/sampling_freq;'])
        
    
        eval(['switch_var = block', num2str(RF_blocks_ind(j)), '.pulse_type;'])
        switch switch_var   %de graaf in vivo spectro
            case 'block'
                eval(['block', num2str(RF_blocks_ind(j)), '.PulseAmpl = block', num2str(RF_blocks_ind(j)), '.flipangle/(settings.general.gamma*block', num2str(RF_blocks_ind(j)), '.dur*settings.TX.block.pulse_integral);'])
                eval(['PulseAmpl_b', num2str(RF_blocks_ind(j)), '_t = block',num2str(RF_blocks_ind(j)), '.PulseAmpl*ones(1, length(timeRF_block', num2str(RF_blocks_ind(j)), '));'])
            case 'sinc'
                %Five lobe sinc (n=3) 
                eval(['B1_norm = sinc(2*3*(timeRF_block', num2str(RF_blocks_ind(j)), '- block', num2str(RF_blocks_ind(j)), '.dur/2)/block', num2str(RF_blocks_ind(j)), '.dur);'])
                eval(['block', num2str(RF_blocks_ind(j)), '.PulseAmpl = block', num2str(RF_blocks_ind(j)), '.flipangle/(settings.general.gamma*block', num2str(RF_blocks_ind(j)), '.dur*settings.TX.sinc.pulse_integral);'])
                eval(['PulseAmpl_b', num2str(RF_blocks_ind(j)), '_t = block', num2str(RF_blocks_ind(j)), '.PulseAmpl*B1_norm;'])
            case 'gaussian10'
                eval(['B1_norm = exp(- 0.1 *(2*(timeRF_block', num2str(RF_blocks_ind(j)), '- block', num2str(RF_blocks_ind(j)), '.dur/2)/block', num2str(RF_blocks_ind(j)), '.dur).^2);'])
                eval(['block', num2str(RF_blocks_ind(j)), '.PulseAmpl = block', num2str(RF_blocks_ind(j)), '.flipangle/(settings.general.gamma*block', num2str(RF_blocks_ind(j)), '.dur*settings.TX.gaussian10.pulse_integral);'])
                eval(['PulseAmpl_b', num2str(RF_blocks_ind(j)), '_t = block', num2str(RF_blocks_ind(j)), '.PulseAmpl*B1_norm;'])
        end
        
        eval(['B1_vec_b', num2str(RF_blocks_ind(j)), '_t = zeros(length(timeRF_block', num2str(RF_blocks_ind(j)), '), matrixsize_signal, matrixsize_signal, 3);'])
        
        for uu = 1:matrixsize_signal
            for uuu = 1:matrixsize_signal
                for ll = 1:3
                    eval(['B1_vec_b', num2str(RF_blocks_ind(j)), '_t(:,',num2str(uu),',', num2str(uuu), ',', num2str(ll),') = 2*PulseAmpl_b', num2str(RF_blocks_ind(j)),'_t.*cos(w*timeRF_block', num2str(RF_blocks_ind(j)),' + block', num2str(RF_blocks_ind(j)),'.initPhase)*B1TX_spatial(', num2str(ll),',', num2str(uu),',', num2str(uuu), ');'])
                end
            end
        end

        eval(['B1_vec_b', num2str(RF_blocks_ind(j)), '_t_st = reshape(B1_vec_b', num2str(RF_blocks_ind(j)), '_t, size(B1_vec_b', num2str(RF_blocks_ind(j)), '_t,1), [], 3);'])
    end

    B1vec_out = cell(length(RF_blocks_ind),1);

    for j = 1:length(RF_blocks_ind)
        eval(['B1vec_out{', num2str(j), '} = B1_vec_b', num2str(RF_blocks_ind(j)), '_t_st;'])
    end
end