function [settings] = InterpolateBFields(settings, ImportB, flag)
%InterpolateBFields Imported magnetic fields are interpolated to desired matrix size
%   Input:  -settings struct
%           -ImportB:  Magnetic field to be interpolated
%           -flag: "Type" of magnetic field: different handling for B0/TX,RX-fields
%   Output:
%           -settings struct

    Nx_B=size(ImportB,2); Ny_B=size(ImportB,3);
    x_B = ([0:Nx_B-1]-Nx_B/2+0.5);y_B = ([0:Ny_B-1]-Ny_B/2+0.5);
    X_B=linspace(x_B(1), x_B(Nx_B), settings.signal.matrixsize_signal);
    Y_B=linspace(y_B(1), y_B(Ny_B), settings.signal.matrixsize_signal);

    X_B_Reco=linspace(x_B(1), x_B(Nx_B), settings.reco.matrixsize_reco);
    Y_B_Reco=linspace(y_B(1), y_B(Ny_B), settings.reco.matrixsize_reco);

    switch flag
        case 'B0'

            BMapInterp = zeros(3,settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);
            BMapInterp_Reco = zeros(3,settings.reco.matrixsize_reco, settings.reco.matrixsize_reco);
            for ll = 1:3
                BMapInterp(ll,:,:) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:)), Y_B', X_B, 'linear');
                BMapInterp_Reco(ll,:,:) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:)), Y_B_Reco', X_B_Reco,'linear');
            end

            settings.signal.B0Map = BMapInterp;
            settings.reco.B0Map = BMapInterp_Reco;
        case 'B1RX'
            settings.RX.nr_receiver = size(ImportB,4);
            BMapInterp = zeros(3,settings.signal.matrixsize_signal, settings.signal.matrixsize_signal, settings.RX.nr_receiver);
            BMapInterp_Reco = zeros(3,settings.reco.matrixsize_reco, settings.reco.matrixsize_reco, settings.RX.nr_receiver);
            for ll = 1:3
                for jj = 1:settings.RX.nr_receiver
                    BMapInterp(ll,:,:,jj) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:,jj)), Y_B', X_B, 'linear');
                    BMapInterp_Reco(ll,:,:,jj) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:,jj)), Y_B_Reco', X_B_Reco, 'linear');
                end
            end

            settings.signal.B1RXMap = BMapInterp;
            settings.reco.B1RXMap = BMapInterp_Reco;
        case 'B1TX'
            BMapInterp = zeros(3,settings.signal.matrixsize_signal, settings.signal.matrixsize_signal);
            BMapInterp_Reco = zeros(3,settings.reco.matrixsize_reco, settings.reco.matrixsize_reco);
            for ll = 1:3
                BMapInterp(ll,:,:) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:)), Y_B', X_B, 'linear');
                BMapInterp_Reco(ll,:,:) =interp2(y_B, x_B, squeeze(ImportB(ll,:,:)), Y_B_Reco', X_B_Reco, 'linear');
            end
            
            settings.signal.B1Map = BMapInterp;
            settings.signal.B1Map_Reco = BMapInterp_Reco;
        otherwise
            warning('Not implemented');
    end
end

