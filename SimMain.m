%% V1: Block wise simulation B-vector 2D
close all; close all hidden; clear all;
addpath(genpath(pwd));                                      %add current path and all subfolders to search path
curr_path = pwd;                                            %needs to be where SimMain is located
settings.general.curr_path = curr_path;

%% definition of pulse sequence to simulate
num_blocks = 3;     %number of total blocks

%block 1
block1.RF_block = true;
block1.readout = false;
if block1.RF_block
    block1.flipangle = 90* pi /180;
    block1.dur = 1*10^-6;%s
    block1.pulse_type = 'block';            %block, sinc, gaussian10

    block1.initPhase = 0;                   %phase with respect to x-axis
    block1.Bfield = '0';
else
    block1.dur = 1*10^-3;%s
    block1.Bfield = '0';
end

%block 2
block2.RF_block = false;
block2.readout = false;
if block2.RF_block
    block2.flipangle = 90* pi /180;
    block2.dur = 100*10^-9;%s
    block2.pulse_type = 'block';            %block, sinc, gaussian10

    block2.initPhase = 0;                   %phase with respect to x-axis
else
    block2.dur = 0.5*10^-3;%s
    block2.Bfield = '1';
end

%block 3
block3.RF_block = false;
block3.readout = true;

if block3.RF_block
    block3.flipangle = 90* pi /180;
    block3.dur = 200*10^-9;%s
    block3.pulse_type = 'block';            %block, sinc, gaussian10

    block3.initPhase = 0;                   %phase with respect to x-axis
else
    block3.dur = 1*10^-3;%s
    block3.Bfield = '2';
end

settings.general.num_blocks = num_blocks;
blocks =cell(num_blocks,1);
for j = 1:num_blocks
    eval(['blocks{', num2str(j), '} = block', num2str(j), ';'])
end
%% parameters
settings.general.B0 = 50*10^-3;%T
settings.reco.FOV =0.2;
settings.signal.matrixsize_signal = 8;
settings.reco.matrixsize_reco = 8;

settings.signal.factor_os = 216;    %Bloch Simulation: Oversampling <--> numerical stability but should be integer of required OS -> easier for downsampling

%factor_os has to be chosen such that factor_os * Nyquist_downsampled = Nyquist_B0Max:
%B0max = settings.general.B0+0.2*settings.trajectory.GRstrength_est
%settings.general.gamma*B0max / (2*pi) / settings.trajectory.BW

settings.trajectory.GRstrength_est = 20*10^-3;
settings.trajectory.type = 'Cart2D';
%% options for simulation
settings.general.relaxation = false;
settings.general.Bloch = false;                                 %otherwise simulated via Matrix-propagation
settings.general.noise = false;                                 %if false: no noise added after signal simulation
settings.general.SNR_dB =  25;

settings.reco.PCG = true;                                      %alternative: ART

settings.general.importB0 = false;
settings.general.importB1_TX = false;
settings.general.importB1_RX = false;
settings.RX.nr_receiver = 1;
settings.general.importEncoding = false;

settings.CST.Nx_CST = 8;                                        %matrix size of field data exported with CST
settings.CST.Ny_CST = settings.CST.Nx_CST;
settings.CST.Nz_CST = settings.CST.Nx_CST;
settings.CST.SixRows = 1;                                       %if 1: data is saved with real and imag parts seperately  -> 6 export fields

%% fundamental constants
settings.general.gamma =  2.675221874411*10^8;            % in rad/s/T; gyromagnetic ratio of H-nuclei.
settings.general.u0=1.26e-6;                              %u_0.
settings.general.T = 293.15;                              %Temperature in K: 20°C
settings.general.k = 1.380649e-23;                        %Boltzmann constant
settings.general.hbar = 1.05457182e-34;                   %hbar, reduced Planck's constant

%RF Pulse LookUp Table
settings.TX.block.pulse_integral = 1;
settings.TX.sinc.pulse_integral = 0.1775;
settings.TX.gaussian10.pulse_integral = 0.565;

[settings] = calcNyquist2D(settings);

%calculation of necessary sampling frequency without mixing the signal
nec_samp_freq = 8*settings.general.gamma*(settings.general.B0 + settings.reco.FOV*settings.trajectory.GRstrength_est)/(2*pi);

disp('Warning: Durations of blocks are overwritten')
%for GRE example: new block lengths definition
scaling_factor = 5;
block2.dur = 1/scaling_factor*settings.trajectory.Tread;
block3.dur = settings.trajectory.Tread;
blocks =cell(num_blocks,1);
for j = 1:num_blocks
    eval(['blocks{', num2str(j), '} = block', num2str(j), ';'])
end

%% coordinate system
FOV = settings.reco.FOV;
matrixsize_signal = settings.signal.matrixsize_signal;
matrixsize_reco = settings.reco.matrixsize_reco;
settings.signal.pixel_width = FOV / matrixsize_signal;
settings.reco.pixel_width = FOV / matrixsize_reco;

%Coordinate System
x = (linspace(-FOV/2, FOV/2- settings.signal.pixel_width, matrixsize_signal));
y = fliplr(x);

[coord_x, coord_y] = meshgrid(x,y);

x_reco = (linspace(-FOV/2, FOV/2- settings.reco.pixel_width, matrixsize_reco));
y_reco = fliplr(x_reco);

[coord_x_reco, coord_y_reco] = meshgrid(x_reco,y_reco);

%% sample
settings.sample.type = 'SheppLogan'; %alternatives: SheppLogan, Sphere, Cylinder, Import, Rect1D
switch settings.sample.type
    case 'SheppLogan'
        
    case 'Sphere'
        settings.sample.sph_shiftx = 2.5*10^-3;         %m
        settings.sample.sph_shifty = -8.52*10^-3;       %m
        settings.sample.sph_shiftz = -8.52*10^-3;       %m
        settings.sample.sph_radius = 0.0725;            %m
        settings.sample.sph_T1 = 0.218;                 %s
        settings.sample.sph_T2 = 0.1735;                %s
    case 'Cylinder'

        settings.sample.cyl_shiftx = -0.05;             %m
        settings.sample.cyl_shifty = 0.01;              %m
        settings.sample.cyl_shiftz = 0.01;              %m
        settings.sample.cyl_radius = 0.04;              %m
        settings.sample.cyl_height = 0.16;              %m
        settings.sample.cyl_T1 = 2.75;                  %s
        settings.sample.cyl_T2 = 2.05;                  %s
    case 'Import'
        settings.sample.path = [curr_path, '/Sample/']; %path where to find the files of a sample
    case 'Rect1D'

    otherwise
        warning('Not implemented');
end

[sampleS, sample_straight] = samplePrep(settings,x,y, coord_x, coord_y);

%% calculation of timings/ time
time_tot_dur =0;
dur_arr = zeros(num_blocks,1);
for j =1:num_blocks
    time_tot_dur = time_tot_dur + blocks{j}.dur;    
    dur_arr(j) = blocks{j}.dur;                     
end
dur_cs = cumsum(dur_arr);

sampling_freq = settings.trajectory.BW*settings.signal.factor_os;

time_tot = 1/sampling_freq:1/sampling_freq:(time_tot_dur+1/sampling_freq);
settings.signal.sampling_freq = sampling_freq;

%indices for individual blocks
indices_start = ones(num_blocks,1);
for j = 2:num_blocks
    indices_start(j) = ceil(dur_cs(j-1)*sampling_freq)+1;
end
indices_end = indices_start(2:end)-1;
indices_end(num_blocks) = length(time_tot);

%% All about Magnetic fields

%B0-Map

B0Map = zeros(3,matrixsize_signal, matrixsize_signal);
%B0Map(3,:,:) =settings.general.B0; 
B0Map_reco = zeros(3,matrixsize_reco,matrixsize_reco);
%B0Map_reco(3,:,:) =settings.general.B0; 

%add imported deltaB0 and add to "Ground"-field
if settings.general.importB0
    pathB0 = [curr_path, '/FieldData/B0Map.mat'];     %path of B0Map: dimensions (3,128,128,128)
    ImportB0 = load(pathB0);
    ImportB0 = ImportB0.B0Map;

    %interpolate to requested size
    [settings] = InterpolateBFields(settings, ImportB0, 'B0');

    B0Map = B0Map + settings.signal.B0Map;
    B0Map_reco = B0Map_reco + settings.reco.B0Map;
    clearvars ImportB0
else

    Desangles = pi/180 * 0.5*repmat((linspace(-40,35, settings.signal.matrixsize_signal)).', [1, settings.signal.matrixsize_signal]);
    
    B0Map(1,:,:) = -settings.general.B0*sin(Desangles);
    B0Map(3,:,:) = settings.general.B0*cos(Desangles);
   
    %Desangles_reco = Desangles;
    Desangles_reco = pi/180 *0.5* repmat((linspace(-40,35, settings.reco.matrixsize_reco)).', [1, settings.reco.matrixsize_reco]);
    
    B0Map_reco(1,:,:) = -settings.general.B0*sin(Desangles_reco);
    B0Map_reco(3,:,:) = settings.general.B0*cos(Desangles_reco);

    anglesB0 = asin(B0Map(1,:,:) / settings.general.B0 )*180/pi;

    %plot field
    % figure;
    % quiver3(x.', ones(matrixsize_signal,1), ones(matrixsize_signal,1), squeeze(B0Map(1,:)).', squeeze(B0Map(2,:)).', squeeze(B0Map(3,:)).');axis equal
    % figure;quiver(x(1),0,squeeze(B0Map(1,1))/1,squeeze(B0Map(3,1)/1));
    % hold on;
    % for j=2:length(B0Map)
    %     quiver(x(j),0,squeeze(B0Map(1,j)/1),squeeze(B0Map(3,j)/1));
    % end
    % hold off;
    %end plot field
end

%B0map_t = permute(repmat(B0map, [1 1 length(time_tot)]), [3,1,2]);

%Transmit B1: Oscillation frequency 'homogeneous' direction may vary at every point in space

if settings.general.importB1_TX
    pathB1 = [curr_path, '/FieldData/B1Map.mat'];     %path of B1Map size: dimensions (3,matrixsize_signal, matrixsize_signal, matrixsize_signal)
    ImportB1 = load(pathB1);
    ImportB1 = ImportB1.B1TX_spatial;

    %interpolate to requested size
    [settings] = InterpolateBFields(settings, ImportB1, 'B1TX');
    B1TX_spatial = settings.signal.B1Map;
    B1TX_spatial_reco = settings.signal.B1Map_Reco;
    clearvars ImportB1
else    
    %if not imported: "build a simulation B1Map"
    B1TX_spatial = zeros(3,matrixsize_signal, matrixsize_signal);
    B1TX_spatial(1,:,:) = 1;

    B1TX_spatial_reco = zeros(3,matrixsize_reco, matrixsize_reco);
    B1TX_spatial_reco(1,:,:) = 1;
end

%time dependence
[settings, B1vec_out, RF_blocks_ind] = TimeDependenceB12D(settings, blocks, indices_start, indices_end, B1TX_spatial, time_tot, sampling_freq);


%Receive B1
if settings.general.importB1_RX
    pathB1RX = [curr_path, '/FieldData/B1RXMap.mat'];     %path of B1RXMap dimensions (3,X,X,X, nr_receiver)
    ImportB1RX = load(pathB1RX);
    ImportB1RX = ImportB1RX.B_RX;

    %interpolate to requested size
    [settings] = InterpolateBFields(settings, ImportB1RX, 'B1RX');
    B_RX = settings.signal.B1RXMap;         %dimensions (3, matrixsize_signal, matrixsize_signal, nr_receiver)
    clearvars ImportB1RX
else
    B_RX = zeros(3, matrixsize_signal,matrixsize_signal, settings.RX.nr_receiver);
    B_RX(2,:,:,:) = 1;
end
B_RX = reshape(B_RX, 3,[]);

%Encoding
Enc_cells = cell(settings.trajectory.N_PhaseEnc, settings.general.num_blocks, 2);    %number of encoding steps, number of sequence blocks, signal / reco

if settings.general.importEncoding
    pathEnc = [curr_path, '/FieldData/Encoding'];
    file_info = natsortfiles(dir(pathEnc));
    filepaths_EncFiles = cell(settings.trajectory.N_PhaseEnc, settings.general.num_blocks);

    startRow = 3;% necessary if there is a header line
    for oo=1:settings.trajectory.N_PhaseEnc
        for ooo=1:settings.general.num_blocks
            filepaths_EncFiles{oo, ooo} = [pathEnc, '/E', num2str(oo), '/B', num2str(ooo), '/E', num2str(oo), '_B', num2str(ooo), '_GRE.txt'];
            [Bout, Bout_Reco] = ImportDataTXT(filepaths_EncFiles{oo,ooo}, settings.CST.SixRows, startRow, settings.CST.Nx_CST, settings.CST.Ny_CST, settings.CST.Nz_CST, settings.signal.matrixsize_signal, settings.reco.matrixsize_reco);
            Enc_cells{oo, ooo,1} = Bout;    %Bout(1:3,X,Y,Z); 1:X-Component, 2:Y-Component, 3:Z-Components
            Enc_cells{oo, ooo,2} = Bout_Reco;
        end
    end
else

    %2D GR cartesian
    for uuu = 1:settings.general.num_blocks
        switch blocks{uuu}.Bfield
            case '0'
                [Enc_cells{:,uuu,1}] = deal(zeros(3, matrixsize_signal, matrixsize_signal));
                [Enc_cells{:,uuu,2}] = deal(zeros(3, matrixsize_reco, matrixsize_reco));
            case '1'%prephase
                amplitudes_y_PE_k = linspace(-settings.trajectory.k_FOV/2, settings.trajectory.k_FOV/2 - settings.trajectory.dK, settings.trajectory.N_PhaseEnc);
                amplitudes_y_PE_G = amplitudes_y_PE_k *2*pi/ (settings.general.gamma * settings.trajectory.Tread);

                for jj = 1:settings.trajectory.N_PhaseEnc
                    B_1_pre = zeros(3, matrixsize_signal, matrixsize_signal);
                    B_1_pre_reco = zeros(3, matrixsize_reco, matrixsize_reco);
                    B_1_pre(3,:,:) = -0.5*scaling_factor*settings.trajectory.GRstrength_est*coord_x + scaling_factor*amplitudes_y_PE_G(jj)*coord_y;
                    %B_1_pre(2,:,:) =(-0.5)*(-0.5)*scaling_factor*settings.trajectory.GRstrength_est*coord_y+ scaling_factor*amplitudes_y_PE_G(jj)*coord_x;     %concomitants
                    B_1_pre_reco(3,:,:) = -0.5*scaling_factor*settings.trajectory.GRstrength_est*coord_x_reco + scaling_factor*amplitudes_y_PE_G(jj)*coord_y_reco;
                    %B_1_pre_reco(2,:,:) = (-0.5)*(-0.5)*scaling_factor*settings.trajectory.GRstrength_est*coord_y_reco + scaling_factor*amplitudes_y_PE_G(jj)*coord_x_reco;
                    Enc_cells{jj,uuu,1} = B_1_pre;
                    Enc_cells{jj,uuu,2} = B_1_pre_reco;

                end
            case '2'%readout
                B_2_pre = zeros(3, matrixsize_signal, matrixsize_signal);
                B_2_pre_reco = zeros(3, matrixsize_reco, matrixsize_reco);
                B_2_pre(3,:,:) = settings.trajectory.GRstrength_est*coord_x;
                %B_2_pre(2,:,:) = -0.5*settings.trajectory.GRstrength_est*coord_y;
                B_2_pre_reco(3,:,:) = settings.trajectory.GRstrength_est*coord_x_reco;
                %B_2_pre_reco(2,:,:) = -0.5*settings.trajectory.GRstrength_est*coord_y_reco;
                [Enc_cells{:,uuu,1}] = deal(B_2_pre);
                [Enc_cells{:,uuu,2}] = deal(B_2_pre_reco);
            otherwise
                error('Option not implemented!')
        end
    end
    %end 2D GR  
end

%% blockwise simulation
tic
Minit = zeros(settings.signal.matrixsize_signal,settings.signal.matrixsize_signal, 3);

for j =1:settings.signal.matrixsize_signal
    for jj = 1:settings.signal.matrixsize_signal
        Minit(j,jj,:) =  sampleS.M0(j,jj)*B0Map(:,j,jj)/norm(B0Map(:,j,jj));    
    end
end


countRO = 0;
for jj = 1:length(blocks)
    if blocks{jj}.readout
        countRO = countRO +1;
        indRO(countRO) = jj;
    end
end

Enc_cells_resh = cell(settings.trajectory.N_PhaseEnc, settings.general.num_blocks,2);
B0Map_st = reshape(B0Map, 3, [],1);
B0Map_reco_st = reshape(B0Map_reco, 3, [],1);

Sig_demod_block = cell(settings.trajectory.N_PhaseEnc, countRO);
Sig_demod_ds_block = cell(settings.trajectory.N_PhaseEnc, countRO);
Sig_demod_ds_EZY_block = cell(settings.trajectory.N_PhaseEnc, countRO);

for jj=1:settings.trajectory.N_PhaseEnc
    countRF = 0;
    counterRO = 0;
    for jjj =1:settings.general.num_blocks
        %determine the necessary magnetic fields
        Enc_cells_resh{jj,jjj,1} = reshape(Enc_cells{jj,jjj,1}, 3, [],1);
        Enc_cells_resh{jj,jjj,2} = reshape(Enc_cells{jj,jjj,2}, 3, [],1);

        B_t_block = permute(repmat((Enc_cells_resh{jj,jjj,1} + B0Map_st).', [1 1 (indices_end(jjj)-indices_start(jjj)+ 1)]), [3,1,2]);
        B_t_block_woRF = permute(repmat((Enc_cells_resh{jj,jjj,1} + B0Map_st).', [1 1 (indices_end(jjj)-indices_start(jjj)+ 1)]), [3,1,2]);
        if blocks{jjj}.RF_block
            countRF = countRF+1;
            B_t_block = B_t_block + B1vec_out{1,countRF};
        end

        %calculate M
        if jjj==1
            M_prev = reshape(Minit, [], 3);
        end

        if settings.general.Bloch
            if settings.general.relaxation
                M = BlochSimMagn_Vector_3D_2DSimRelax_block(settings, time_tot(indices_start(jjj):indices_end(jjj)), B_t_block, sampleS, M_prev);
            else
                M = BlochSimMagn_Vector_3D_2DSim_block(settings, time_tot(indices_start(jjj):indices_end(jjj)), B_t_block, M_prev);
            end
        else
            if blocks{jjj}.RF_block
                M = ExcMatrixCalc_block(settings, B_t_block_woRF, indices_start, indices_end, B1TX_spatial, blocks, jjj, M_prev, false);
            else
                if settings.general.relaxation
                    M = PropMatrixBconstRelax_block(settings, indices_start, indices_end, B_t_block_woRF, time_tot, M_prev, sampleS, Minit, jjj);    
                else
                    M = PropMatrixBconstSum_block(settings, indices_start, indices_end, B_t_block_woRF, 1, M_prev, jjj);

                end
            end
        end

        M_prev = squeeze(M(end,:,:));

        %calculate Signal from M
        if blocks{jjj}.readout
            counterRO = counterRO+1;

            %downsample M
            M = squeeze(M(1:round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW)):end,:,:));
            time_RO_ds = squeeze(time_tot(indices_start(jjj):round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW)):indices_end(jjj)));
            V= zeros(1,length(time_RO_ds));

            for ll = 1:length(time_RO_ds)
                V(ll) = -1*sum(squeeze(B_RX(1,:).').*squeeze(M(ll, :,1)).'+ squeeze(B_RX(2,:).').*squeeze(M(ll, :,2)).' + squeeze(B_RX(3,:).').*squeeze(M(ll, :,3)).');
            end
            V = gradient(squeeze(V)) ./ gradient(time_RO_ds);

            demod_freq = (settings.general.gamma*settings.general.B0);
            demod_sig = exp(-1i*demod_freq*time_RO_ds);

            sig_demod_NPE = hilbert(V).*demod_sig;
            sig_demod_NPE_dec = decimate(sig_demod_NPE, round(nec_samp_freq/settings.trajectory.BW));
            sig_Demod_NPE_decEZY = squeeze(sig_demod_NPE(2:round(nec_samp_freq/settings.trajectory.BW):end));

            Sig_demod_block{jj,counterRO} = sig_demod_NPE;
            Sig_demod_ds_block{jj,counterRO} = sig_demod_NPE_dec;
            Sig_demod_ds_EZY_block{jj,counterRO} = sig_Demod_NPE_decEZY;
        end
    end
end

time_tot_ds = time_tot(1:settings.signal.factor_os:end);
time_tot = time_tot(1:round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW)):end);
sampling_freq_M = nec_samp_freq;
sampling_freq_ds = settings.trajectory.BW;

%reorder signal cells
Sig_demod_cell = cell(countRO,1);
Sig_demod_ds_cell = cell(countRO,1);
Sig_demod_dsEZY_cell = cell(countRO,1);

for jj = 1:countRO
    mat_temp = cell2mat(Sig_demod_block(:,jj));
    Sig_demod_cell{jj} = mat_temp;

    mat_temp = cell2mat(Sig_demod_ds_block(:,jj));
    Sig_demod_ds_cell{jj} = mat_temp;

    mat_temp = cell2mat(Sig_demod_ds_EZY_block(:,jj));
    Sig_demod_dsEZY_cell{jj} = mat_temp;
end

if settings.general.noise
    for j = 1:countRO
        Sig_demod_cell{j} =  awgn((Sig_demod_cell{j}), settings.general.SNR_dB, 'measured');    %add Noise such that we achieve the desired SNR
        Sig_demod_ds_cell{j} =  awgn((Sig_demod_ds_cell{j}), settings.general.SNR_dB, 'measured'); 
        Sig_demod_dsEZY_cell{j} = awgn((Sig_demod_dsEZY_cell{j}), settings.general.SNR_dB, 'measured'); 
    end
end
toc
%clearvars M V Sig_demod_ds_cell
clearvars amplitudes_y_PE_G amplitudes_y_PE_k anglesB0 B0Map B0Map_reco B0Map_st B1TX_spatial B1vec_out B_1_pre B_1_pre_reco B_2_pre B_2_pre_reco B_t_block B_t_block_woRF coord_x coord_x_reco coord_y coord_y_reco demod_sig Desangles Desangles_reco dur_cs dur_arr Enc_cells mat_temp Minit Sig_demod_block Sig_demod_ds_block Sig_demod_ds_EZY_block sig_demod_NPE sig_demod_NPE_dec sig_Demod_NPE_decEZY time_RO_ds x x_reco y y_reco sample_straight

%% Reconstruction
tic
%calculate M at end of each block
M_endBlock = cell(settings.general.num_blocks,1);

e_B0_reco_init = zeros(settings.reco.matrixsize_reco^2, 3, settings.trajectory.N_PhaseEnc);

for j =1:settings.reco.matrixsize_reco^2
    for jjj = 1:settings.trajectory.N_PhaseEnc
        e_B0_reco_init(j,:,jjj) = B0Map_reco_st(:,j)/norm(B0Map_reco_st(:,j));
    end
end

%calculate magnetic fields
B_t_woRF_reco_block = zeros(settings.reco.matrixsize_reco^2, 3,settings.trajectory.N_PhaseEnc,settings.general.num_blocks);
for j = 1:settings.general.num_blocks
    for jj=1:settings.trajectory.N_PhaseEnc
        B_t_woRF_reco_block(:,:,jj,j) = (Enc_cells_resh{jj,j,2} + B0Map_reco_st).'; 
    end
end
clearvars B0Map_reco_st Enc_cells_resh
for j = 1:settings.general.num_blocks
    
    if blocks{j}.RF_block
        if j ==1
            M_prop_reco = ExcMatrixCalc_Reco(settings, squeeze(B_t_woRF_reco_block(:,:,:,j)), indices_start, indices_end, B1TX_spatial_reco, blocks, j, e_B0_reco_init, true);
            M_endBlock{j} = squeeze(M_prop_reco(end,:,:,:));
        else
            M_prev = M_endBlock{j-1};
            M_prop_reco = ExcMatrixCalc_Reco(settings, squeeze(B_t_woRF_reco_block(:,:,:,j)), indices_start, indices_end, B1TX_spatial_reco, blocks, j, M_prev, true);
            M_endBlock{j} = squeeze(M_prop_reco(end,:,:,:));
        end
    else
        M_prop_reco = zeros(settings.reco.matrixsize_reco^2,3,settings.trajectory.N_PhaseEnc);
        for uuu =1:settings.trajectory.N_PhaseEnc
            for ll = 1:settings.reco.matrixsize_reco^2
                B_block_reco = squeeze(B_t_woRF_reco_block(ll,:,uuu,j));
                u_B_block_reco = squeeze(B_block_reco) / norm(squeeze(B_block_reco));%unit vector during different blocks
            
                x1_block_r = 1/sqrt(2*(u_B_block_reco(1)^2 + u_B_block_reco(3)^2))*[u_B_block_reco(1)*u_B_block_reco(2) - 1i*u_B_block_reco(3) ; - (u_B_block_reco(1)^2 + u_B_block_reco(3)^2) ; u_B_block_reco(2)*u_B_block_reco(3) + 1i*u_B_block_reco(1)];
                x2_block_r = u_B_block_reco.';
                x3_block_r = conj(x1_block_r);
    
                B_mag_reco = norm(B_block_reco);
                lambda3_r_t = 1i*settings.general.gamma*(B_mag_reco - settings.general.B0)*blocks{j}.dur;
                lambda1_r_t = conj(lambda3_r_t);
                lambda2_r_t = 0;
    
                if j ==1
                    M_prev = e_B0_reco_init;
                else
                    M_prev = M_endBlock{j-1};
                end
                M_prop_reco(ll,:,uuu)= exp(lambda1_r_t)*x1_block_r'*M_prev(ll,:,uuu).'*x1_block_r + exp(lambda2_r_t)*x2_block_r'*M_prev(ll,:,uuu).'*x2_block_r + exp(lambda3_r_t)*x3_block_r'*M_prev(ll,:,uuu).'*x3_block_r;
            end
        end
        M_endBlock{j} = M_prop_reco;
    end
end
clearvars M_prop_reco M_prev lambda1_r_t lambda2_r_t lambda3_r_t B_mag_reco
if settings.reco.PCG
%calculations of relevant contributions during RO
disp('TO BE CHANGED with multiple ROs: indRO');
Enc_r = cell(countRO,settings.trajectory.N_PhaseEnc);

for j = 1:countRO
    for uuu=1:settings.trajectory.N_PhaseEnc
    rx_ev_r = zeros(matrixsize_reco^2,1);
    ev_imagn_r = zeros(matrixsize_reco^2,1);
    B_tp_mag_r_t = zeros(matrixsize_reco^2,1);
    lambda3_r_t = zeros((indices_end(indRO(j))-indices_start(indRO(j))+1)/round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW)), matrixsize_reco^2);
    for ll =1:matrixsize_reco^2

        B_RO_reco = squeeze(B_t_woRF_reco_block(ll,:,uuu, indRO)); 
        u_B_RO_reco = squeeze(B_RO_reco) / norm(squeeze(B_RO_reco));%unit vector during Readout
        
        x1_r = 1/sqrt(2*(u_B_RO_reco(1)^2 + u_B_RO_reco(3)^2))*[u_B_RO_reco(1)*u_B_RO_reco(2) - 1i*u_B_RO_reco(3) ; - (u_B_RO_reco(1)^2 + u_B_RO_reco(3)^2) ; u_B_RO_reco(2)*u_B_RO_reco(3) + 1i*u_B_RO_reco(1)];
        x2_r = u_B_RO_reco;
        x3_r = conj(x1_r);

        rx_ev_r(ll) = squeeze(B_RX(:,ll))'*x3_r;
        if indRO(j) ==1
            M_prev_reco_RO = e_B0_reco_init;
        else
            M_prev_reco_RO = M_endBlock{indRO(j)-1};
        end
        ev_imagn_r(ll) = x3_r'*squeeze(M_prev_reco_RO(ll,:,uuu).');

        B_tp_mag_r_t(ll) = norm(B_RO_reco);
        lambda3_r_t(:,ll) = (1i*settings.general.gamma*(B_tp_mag_r_t(ll)-settings.general.B0)*(1:(indices_end(indRO(j))-indices_start(indRO(j))+1)/(round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW))))*1/sampling_freq_M).';
    end
             
    Enc_r{j,uuu} = -1i*settings.general.gamma*(B_tp_mag_r_t.*rx_ev_r.*ev_imagn_r).'.*exp(lambda3_r_t);
    end
end
clearvars e_B0_reco_init rx_ev_r ev_imagn_r B_tp_mag_r_t lambda3_r_t M_prev_reco_RO

E_r = cat(1, Enc_r{1,:});
signal_r = Sig_demod_cell{1}.';
signal_r = signal_r(:);

for jj = 1:settings.trajectory.N_PhaseEnc
    for ll = 1:matrixsize_reco^2
        E_PE_indv = Enc_r{1,jj};
        E_PE_indv_dec(:,ll)= decimate(squeeze(E_PE_indv(:,ll)), round(nec_samp_freq/settings.trajectory.BW));
        E_PE_indv_decEZY(:,ll)= squeeze(E_PE_indv(2:round(nec_samp_freq/settings.trajectory.BW):end,ll));
    end
    E_PE_r{1,jj} = E_PE_indv_decEZY; 
end
E_r_RO_dec = cat(1, E_PE_r{1,:});

iteration = 100000;

% wS = 1;
% M_SSOR = wS/(2-wS)*(1/wS*diag(diag(E_r'*E_r)) + tril(E_r'*E_r))*(diag(1/diag(E_r'*E_r)))*(1/wS*diag(diag(E_r'*E_r)) + triu(E_r'*E_r));
reco_rho_straight_rrr = pcg(E_r'*E_r, E_r'*signal_r, 1e-15, iteration, diag(diag(E_r'*E_r)));%+0.001*eye(size(E_r'*E_r))
reco_rho_straight_rrr_img = reshape(reco_rho_straight_rrr, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, []);

reco_rho_ttt= E_r \ signal_r;
reco_rho_ttt_img = reshape(reco_rho_ttt, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, []);

Sig_demod_ds = Sig_demod_dsEZY_cell{1}.';
Sig_demod_ds = Sig_demod_ds(:);

reco_rho_straight_rrr_nyq = pcg(E_r_RO_dec'*E_r_RO_dec, E_r_RO_dec'*Sig_demod_ds, 1e-15, iteration);%*signal_r
reco_rho_straight_rrr_nyq_img = reshape(reco_rho_straight_rrr_nyq, settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, []);

else
%% ART reco

%Parameter that worked well
max_it_ART = 16;
lambda = 0.08; 
N_ART = settings.reco.matrixsize_reco^2;
x0_ART = zeros(N_ART,1);
signal_r = Sig_demod_cell{1}.';
signal_r = signal_r(:);
S_backStraightCut = signal_r;
time_length = (indices_end(indRO(1))-indices_start(indRO(1))+1)/round(settings.signal.factor_os/(nec_samp_freq / settings.trajectory.BW));

if false
    B_t_woRF_reco_block = gpuArray((B_t_woRF_reco_block));
    S_backStraightCut = gpuArray(S_backStraightCut);
    X = gpuArray(x0_ART);
    e_B0_reco_init = gpuArray(e_B0_reco_init);
    M_endBlock_gpu = gpuArray(M_endBlock{2});
    B_RX = gpuArray(B_RX);
    rx_ev_r = gpuArray(zeros(matrixsize_reco^2,1));
    ev_imagn_r = gpuArray(zeros(matrixsize_reco^2,1));
    B_tp_mag_r_t = gpuArray(zeros(matrixsize_reco^2,1));
    lambda3_r_t = gpuArray(zeros(1,matrixsize_reco^2));
else
    X = (x0_ART);
end

for j = 1:max_it_ART
    for jj=1:length(S_backStraightCut)

        if mod(jj,time_length)==0
            t_index = time_length;
        else
            t_index = mod(jj, time_length);
        end

        if mod(jj,time_length)==0
            N_PE_index = floor(jj/time_length);
        else
            N_PE_index = floor(jj/time_length) +1;
        end

        rx_ev_r = (zeros(matrixsize_reco^2,1));
        ev_imagn_r = (zeros(matrixsize_reco^2,1));
        B_tp_mag_r_t = (zeros(matrixsize_reco^2,1));
        lambda3_r_t = (zeros(1,matrixsize_reco^2));

        for ll =1:matrixsize_reco^2
            B_RO_reco = squeeze(B_t_woRF_reco_block(ll,:,N_PE_index, indRO)); 
            u_B_RO_reco = squeeze(B_RO_reco) / norm(squeeze(B_RO_reco));%unit vector during Readout

            x1_r = 1/sqrt(2*(u_B_RO_reco(1)^2 + u_B_RO_reco(3)^2))*[u_B_RO_reco(1)*u_B_RO_reco(2) - 1i*u_B_RO_reco(3) ; - (u_B_RO_reco(1)^2 + u_B_RO_reco(3)^2) ; u_B_RO_reco(2)*u_B_RO_reco(3) + 1i*u_B_RO_reco(1)];
            x2_r = u_B_RO_reco;
            x3_r = conj(x1_r);

            rx_ev_r(ll) = squeeze(B_RX(:,ll))'*x3_r;
            if indRO(1) ==1
                M_prev_reco_RO = e_B0_reco_init;
            else
                M_prev_reco_RO = M_endBlock{indRO(1)-1};
                %M_prev_reco_RO = M_endBlock_gpu;
            end
            ev_imagn_r(ll) = x3_r'*squeeze(M_prev_reco_RO(ll,:,N_PE_index).');

            B_tp_mag_r_t(ll) = norm(B_RO_reco);
            lambda3_r_t(ll) = 1i*settings.general.gamma*(B_tp_mag_r_t(ll) - settings.general.B0)*t_index*1/sampling_freq_M;
        end

        An = -1i*settings.general.gamma*B_tp_mag_r_t.*rx_ev_r.*ev_imagn_r.*exp(lambda3_r_t).';
        An = squeeze(An).';

        unitLen = norm(An);
        d = (S_backStraightCut(jj)-An*X)./unitLen^2;
        X = X+(lambda*d.*conj(An)).';
    end
end

reco_rho_img_ART = reshape(gather(X), settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, []);
end
toc