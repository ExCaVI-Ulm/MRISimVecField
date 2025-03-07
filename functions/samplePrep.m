function [sampleS, sample_straight] = samplePrep(settings,x,y,coord_x, coord_y)
%samplePrep All about sample/phantom definition
%   Input:  -settings struct
%           -x,y,coord_x, coord_y: coordinate system: x,y linear 1D;
%           coord_x, coord_y 2D coordinates
%   Output:
%           -sampleS: sample struct
%           -sample_straight: same as sampleS but in matrix form and as 1D
%           vector

disp('Sample');tic

FOV = settings.reco.FOV;
matrixsize_signal = settings.signal.matrixsize_signal;
matrixsize_reco = settings.reco.matrixsize_reco;

switch settings.sample.type
    case 'SheppLogan'

        shepplogan = (phantom('Modified Shepp-Logan', matrixsize_reco));
        sampleS.M0_reco = shepplogan;
        shepplogan = imresize(shepplogan,[matrixsize_signal matrixsize_signal],'nearest');
        sampleS.M0 = shepplogan;

        sample_straight = reshape(shepplogan, 1, []);
        sample_straight(2,:) =  reshape(coord_x, 1, []); %x
        sample_straight(3,:) =  reshape(coord_y, 1, []); %y

        %T2

        T2d = [  70 .69   .92    0     0     0   
                100-70  .6624 .8740   0  -.0184   0
                500-100 .1600 .4100 -.22    0     18
                500-100  .1100 .3100  .22    0    -18
                120-100  .2100 .2500   0    .35    0
                60-100  .0460 .0460   0    .1     0
                120-100  .0460 .0460   0   -.1     0
                80-100  .0460 .0230 -.08  -.605   0 
                80-100  .0230 .0230   0   -.606   0
                80-100  .0230 .0460  .06  -.605   0   ];

        %T1
        T1d = [  450  .69   .92    0     0     0   
                950-450  .6624 .8740   0  -.0184   0
                2000-950  .1100 .3100  .22    0    -18
                2000-950 .1600 .4100 -.22    0     18
                600-950  .2100 .2500   0    .35    0
                1200-950  .0460 .0460   0    .1     0
                1200-950  .0460 .0460   0   -.1     0
                1200-950  .0460 .0230 -.08  -.605   0 
                1200-950  .0230 .0230   0   -.606   0
                1200-950  .0230 .0460  .06  -.605   0   ];
        %T2
        P2 = phantom(T2d,matrixsize_signal);
        P2(P2==0) = Inf;    %T2 = 0 is replaced by Inf
        sampleT2 = P2*10^-3;
        sampleS.T2 = sampleT2;
        sample_straight(5, :) = reshape(sampleT2, 1, []); %T2

        %T1
        P1 = phantom(T1d,matrixsize_signal);
        P1(P1==0) = Inf;    %T1 = 0 is replaced by Inf
        sampleT1 = P1*10^-3;
        sampleS.T1 = sampleT1;
        sample_straight(6, :) = reshape(sampleT1, 1, []); %T1
        
    case 'PointPh'
        size_points = 4;%px
        number_pf_points = matrixsize_reco /size_points;

        subblocks = [0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0;0 0 1 1 1 1 0 0;0 0 1 1 1 1 0 0;0 0 1 1 1 1 0 0;0 0 1 1 1 1 0 0;0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 0];
        conc_d1 = cat(1, subblocks, subblocks, subblocks,subblocks,subblocks,subblocks, subblocks, subblocks);
        M0 = cat(2,conc_d1, conc_d1, conc_d1,conc_d1,conc_d1,conc_d1, conc_d1, conc_d1);
        sampleS.M0_reco = M0;
        sampleS.M0 = imresize(M0,[matrixsize_signal matrixsize_signal],'nearest');

        sample_straight = reshape(M0, 1, []);
        sample_straight(2,:) =  reshape(coord_x, 1, []); %x
        sample_straight(3,:) =  reshape(coord_y, 1, []); %y
        

        sampleS.T1 = settings.sample.cyl_T1*squeeze(M0);
        sampleS.T1(sampleS.T1==0) = Inf;    %T1 = 0 is replaced by Inf --> leads to wrong values for scaling of T1 but there is no Magnetization either
        sampleS.T2 = settings.sample.cyl_T2*squeeze(M0);
        sampleS.T2(sampleS.T2==0) = Inf;    %T2 = 0 is replaced by Inf
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        sampleS.DB = zeros(matrixsize_signal, matrixsize_signal);

    case 'Sphere'
        SpherePhantom = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        
        x_s = x;
        y_s = y;
        z_s = z;
        [coord_s_x, coord_s_y, coord_s_z] = meshgrid(x_s,y_s, z_s);
        
        I_sphPh=(coord_s_x-settings.sample.sph_shiftx).^2+(coord_s_y-settings.sample.sph_shifty).^2 +(coord_s_z-settings.sample.sph_shiftz).^2< settings.sample.sph_radius^2; %
        
        SpherePhantom(I_sphPh)=1;%M0        
        sampleS.M0 = squeeze(SpherePhantom);
        
        sampleS.T1 = settings.sample.sph_T1*squeeze(SpherePhantom);
        sampleS.T1(sampleS.T1==0) = Inf;    %T1 = 0 is replaced by Inf
        sampleS.T2 = settings.sample.sph_T2*squeeze(SpherePhantom);
        sampleS.T2(sampleS.T2==0) = Inf;    %T2 = 0 is replaced by Inf

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(coord_s_x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(coord_s_y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(coord_s_z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
        Susc_sph = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        Susc_sph(I_sphPh) = abs(settings.sample.sph_SusceptInt - settings.sample.sph_SusceptExt);
        
        clearvars coord_s_x coord_s_y coord_s_z Susc_sph SpherePhantom x_s y_s z_s I_sphPh
    case 'Cylinder'
        CylinderPhantom = zeros(matrixsize_signal, matrixsize_signal, matrixsize_signal);
        bool_cyl = (and((coord_z-settings.sample.cyl_shiftz).^2 + (coord_y-settings.sample.cyl_shifty).^2 <= settings.sample.cyl_radius^2, and(coord_x>settings.sample.cyl_shiftx,coord_x<(settings.sample.cyl_height+settings.sample.cyl_shiftx))));  %&& Z>0 && Z<height_cylinder
        
        CylinderPhantom(bool_cyl)=1;%M0        
        sampleS.M0 = squeeze(CylinderPhantom);
        
        sampleS.T1 = settings.sample.cyl_T1*squeeze(CylinderPhantom);
        sampleS.T1(sampleS.T1==0) = Inf;    %T1 = 0 is replaced by Inf --> leads to wrong values for scaling of T1 but there is no Magnetization either
        sampleS.T2 = settings.sample.cyl_T2*squeeze(CylinderPhantom);
        sampleS.T2(sampleS.T2==0) = Inf;    %T2 = 0 is replaced by Inf

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(coord_x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(coord_y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(coord_z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1        

    case 'Import'
        %examplary load JEMRIS Brain Phantom
        load([settings.sample.path,'MNIbrain.mat'], 'BRAIN');
        load([settings.sample.path,'MNIdeltaB.mat'], 'DB');

        %Code from JEMRIS

        %tissue parameters from JEMRIS, mrm.29009 (mean) @50mT, Bottomley NMR Relax in tissue
        %assuming: T2 is independent of field strength
        %        T1  T2 T2*[ms]  M0 CS[rad/sec]      Label
        tissue=[3702 329  158   1.00   0         ;  % 1 = CSF
                 326.7  83   69   0.86   0         ;  % 2 = GM
                 274.7  70   61   0.77   0         ;  % 3 = WM
                 146  70   58   1.00 7.9*2*pi    ;  % 4 = Fat (CS @ 53mT Tesla) ->suscept map
                 130  47   30   1.00   0         ;  % 5 = Muscle / Skin
                 450 329   58   1.00   0         ;  % 6 = Skin
                   0   0    0   0.00   0         ;  % 7 = Skull
                 826  83   69   0.86   0         ;  % 8 = Glial Matter
                 122  70   61   0.77   0         ;];% 9 = Meat

        %parameter maps
        PARAMS={'M0','T1','T2','T2S','DB'};
        fact=[1 1 1 1 1]; %if 1: activated; M0, T1, T2, T2s, CS
        INDEX =[4 1 2 3 5];
        for i=1:9
            for j=1:5
                if i==1,eval(['BrainSample.',PARAMS{j},'=zeros(size(BRAIN));']);end
                I   = find(BRAIN==i);
                ind = INDEX(j);
                eval(['BrainSample.',PARAMS{j},'(I)=fact(j)*tissue(i,ind);']);
            end
        end

        %add susceptibility issues
        %interpolate DB to twice the initial size
        Nx_suscep=size(DB,1); x_suscep=([0:Nx_suscep-1]-Nx_suscep/2+0.5);
        Ny_suscep=size(DB,2); y_suscep=([0:Ny_suscep-1]-Ny_suscep/2+0.5);
        Nz_suscep=size(DB,3); z_suscep=([0:Nz_suscep-1]-Nz_suscep/2+0.5);
        X_suscep=[x_suscep(1):(x_suscep(Nx_suscep)-x_suscep(1))/(2*Nx_suscep-1):x_suscep(Nx_suscep)];
        Y_suscep=[y_suscep(1):(y_suscep(Ny_suscep)-y_suscep(1))/(2*Ny_suscep-1):y_suscep(Ny_suscep)];
        Z_suscep=[z_suscep(1):(z_suscep(Nz_suscep)-z_suscep(1))/(2*Nz_suscep-1):z_suscep(Nz_suscep)];

        DB = interp3(y_suscep,x_suscep,z_suscep,DB,Y_suscep',X_suscep,Z_suscep,'spline');

        BrainSample.DB = BrainSample.DB ;%+ 2*pi*1e6*DB*settings.general.FreqField;
        for j =1:5
            eval(['BrainSample.',PARAMS{j},'=flip(permute(shiftdim(BrainSample.',PARAMS{j},',2),[1,3,2]),1);']);
        end

        %interpolate sample to size of B0
        Nx_interB_sample=size(BrainSample.M0,1); Ny_interB_sample=size(BrainSample.M0,2); Nz_interB_sample=size(BrainSample.M0,3);
        x_interB_sample = ([0:Nx_interB_sample-1]-Nx_interB_sample/2+0.5);y_interB_sample = ([0:Ny_interB_sample-1]-Ny_interB_sample/2+0.5);z_interB_sample = ([0:Nz_interB_sample-1]-Nz_interB_sample/2+0.5);
        X_interB_sample=linspace(x_interB_sample(1), x_interB_sample(Nx_interB_sample), matrixsize_signal);
        Y_interB_sample=linspace(y_interB_sample(1), y_interB_sample(Ny_interB_sample), matrixsize_signal);
        Z_interB_sample=linspace(z_interB_sample(1), z_interB_sample(Nz_interB_sample), matrixsize_signal);
        
        for j =1:5
            eval(['BrainSample.',PARAMS{j},'=interp3(y_interB_sample, x_interB_sample, z_interB_sample, BrainSample.',PARAMS{j},', Y_interB_sample'', X_interB_sample, Z_interB_sample, ''spline'');']);
        end
        
        M03D = BrainSample.M0;
        T13D = BrainSample.T1;
        T23D = BrainSample.T2;
        T2S3D = BrainSample.T2S;
        dB3D = BrainSample.DB;
        BrainSample.M0 = M03D(:,:,ceil(size(M03D,3)/2));
        BrainSample.T1 = T13D(:,:,ceil(size(M03D,3)/2));
        BrainSample.T2 = T23D(:,:,ceil(size(M03D,3)/2));
        BrainSample.T2S = T2S3D(:,:,ceil(size(M03D,3)/2));
        BrainSample.DB = 0*dB3D(:,:,ceil(size(M03D,3)/2));

        sampleS.M0 = squeeze(BrainSample.M0);
        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(coord_x, 1, []); %x
        sample_straight(3,:) =  reshape(coord_y, 1, []); %y
        %sample_straight(4,:) =  reshape(coord_z, 1, []); %z

        %T2
        T2_jemB = BrainSample.T2;
        T2_jemB(T2_jemB<eps) = Inf;    %T2 = 0 is replaced by Inf
        T2_jemB = T2_jemB*10^-3;
        sampleS.T2 = T2_jemB;
        sample_straight(5, :) = reshape(T2_jemB, 1, []); %T2

        %T1
        T1_jemB = BrainSample.T1;
        T1_jemB(T1_jemB<eps) = Inf;    %T1 = 0 is replaced by Inf
        T1_jemB = T1_jemB*10^-3;
        sampleS.T1 = T1_jemB;
        sample_straight(6, :) = reshape(T1_jemB, 1, []); %T1
        

    case 'Rect1D'
        RectPhantom = rectangularPulse(-settings.reco.FOV/3,settings.reco.FOV/3,x);
          
        sampleS.M0 = squeeze(RectPhantom);
        
        sampleS.T1 = 700*10^-3*ones(size(RectPhantom));%Inf(size(RectPhantom));
        sampleS.T2 = 100*10^-3*ones(size(RectPhantom));%Inf(size(RectPhantom));

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(y), 1, []); %y
        sample_straight(4,:) =  reshape(squeeze(z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
    case 'Rect2D'
        RectPhantom = repmat(rectangularPulse(-settings.reco.FOV/3,settings.reco.FOV/3,x), [length(x),1]);
          
        sampleS.M0 = squeeze(RectPhantom);
        
        sampleS.T1 = 700*10^-3*ones(size(RectPhantom));%Inf(size(RectPhantom));
        sampleS.T2 = 100*10^-3*ones(size(RectPhantom));%Inf(size(RectPhantom));

        sample_straight = reshape(sampleS.M0, 1, []);
        sample_straight(2,:) =  reshape(squeeze(coord_x), 1, []); %x
        sample_straight(3,:) =  reshape(squeeze(coord_y), 1, []); %y
        sample_straight(4,:) =  0;%reshape(squeeze(z), 1, []); %z
        sample_straight(5, :) = reshape(sampleS.T2, 1, []); %T2
        sample_straight(6, :) = reshape(sampleS.T1, 1, []); %T1
        
    otherwise
        warning('Not implemented yet!')
end

as(sampleS.M0)
timeElapsed_sample = toc
end