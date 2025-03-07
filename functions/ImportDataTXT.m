function [Bout, Bout_Reco] = ImportDataTXT(filepath, SixRows, startRow, Nx_CST, Ny_CST, Nz_CST, matrixsize_signal, matrixsize_reco)
%ImportDataTXT Import of magnetic fiedl data
%   Input:  -filepath: path where files to be imported are located
%           -SixRows:  Specifies format of data to be read in (3 or 6 magnetic field components)
%           -Nx_CST, Ny_CST, Nz_CST: Size of matrix being imported from CST
%           -matrixsize_signal, matrixsize_reco: matrixsizes either for
%           signal or reconstruction
%   Output:
%           -Bout: Imported magnetic vector field
%           -Bout_Reco: Imported magnetic vector field for reconstruction

    if SixRows
    	    %specify format which should be read in
    
            %formatSpec = '%*51s%17f%17f%17f%17f%f%[^\n\r]';
            %formatSpec = [repmat('%*s',1,3), repmat('%f',1,6)];%skip first 3 columns <-> coordinates, read 6 values
            formatSpec = [repmat('%*s',1,2), repmat('%f',1,3)];
            SixRows = false;
    else
        formatSpec = '%*51s%17f%17f%f%[^\n\r]';
    end

    fileID = fopen(filepath,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', ',', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-3, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    
    fclose(fileID);
    if SixRows
        BXC1 = dataArray{:, 1};
        BYC1 = dataArray{:, 3};
        BZC1 = dataArray{:, 5};
    else
        BXC1 = dataArray{:, 1};
        BYC1 = dataArray{:, 2};
        BZC1 = dataArray{:, 3};
    end

    BXrFOV = flip(reshape(BXC1, Nx_CST,Ny_CST),2).';
    BYrFOV = flip(reshape(BYC1, Nx_CST,Ny_CST),2).';
    BZrFOV = flip(reshape(BZC1, Nx_CST,Ny_CST),2).';

    %interpolate to desired matrixsize
    
    Nx_Bfield=size(BXrFOV,1); Ny_Bfield=size(BXrFOV,2);
    x_Bfield = ([0:Nx_Bfield-1]-Nx_Bfield/2+0.5);y_Bfield = ([0:Ny_Bfield-1]-Ny_Bfield/2+0.5);
    X_Bfield=linspace(x_Bfield(1), x_Bfield(Nx_Bfield), matrixsize_signal);
    Y_Bfield=linspace(y_Bfield(1), y_Bfield(Ny_Bfield), matrixsize_signal);

    X_Bfield_Reco=linspace(x_Bfield(1), x_Bfield(Nx_Bfield), matrixsize_reco);
    Y_Bfield_Reco=linspace(y_Bfield(1), y_Bfield(Ny_Bfield), matrixsize_reco);

    BXrFOV_interpol=interp2(y_Bfield, x_Bfield, BXrFOV,Y_Bfield', X_Bfield, 'linear');
    BYrFOV_interpol=interp2(y_Bfield, x_Bfield, BYrFOV,Y_Bfield', X_Bfield, 'linear');
    BZrFOV_interpol=interp2(y_Bfield, x_Bfield, BZrFOV,Y_Bfield', X_Bfield, 'linear');

    BXrFOV_interpol_Reco=interp2(y_Bfield, x_Bfield, BXrFOV,Y_Bfield_Reco', X_Bfield_Reco, 'linear');
    BYrFOV_interpol_Reco=interp2(y_Bfield, x_Bfield, BYrFOV,Y_Bfield_Reco', X_Bfield_Reco, 'linear');
    BZrFOV_interpol_Reco=interp2(y_Bfield, x_Bfield, BZrFOV,Y_Bfield_Reco', X_Bfield_Reco, 'linear');

    BFOV_interpol_Reco=sqrt(BXrFOV_interpol_Reco.^2 + BYrFOV_interpol_Reco.^2 + BZrFOV_interpol_Reco.^2);

    Bout = zeros(3, matrixsize_signal, matrixsize_signal);
    Bout(1,:,:) = BXrFOV_interpol;    
    Bout(2,:,:) = BYrFOV_interpol;
    Bout(3,:,:) = BZrFOV_interpol;

    Bout_Reco = zeros(3, matrixsize_reco, matrixsize_reco);
    Bout_Reco(1,:,:) = BXrFOV_interpol_Reco;
    Bout_Reco(2,:,:) = BYrFOV_interpol_Reco;
    Bout_Reco(3,:,:) = BZrFOV_interpol_Reco;
end