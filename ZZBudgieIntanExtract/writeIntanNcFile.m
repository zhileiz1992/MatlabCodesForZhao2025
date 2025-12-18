function writeIntanNcFile(filepath, timeVector, deltaT, channel, metaData, data, overwrite)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% writeIntanNcFile: A function to create a binary netCDF format file from
%   a 1 channel time series, for example a single channel of Intan data.
%
% usage:  
%   writeIntanNcFile(filepath, timeVector, deltaT, channel, metaData, data
%                    [, overwrite])
%
% where,
%    filepath is a char array representing the path to save the binary file
%       to. If the filepath lacks a file extension, '.nc' will be added.
%    timeVector is a time/date vector of the form
%           [year, month, day, hour, minute, second, microseconds]
%       or
%           [year, month, day, hour, minute, fractionalSeconds]
%       or
%           A numerical serial datenum (see datenum function)
%    deltaT is a number indicating the sampling period in seconds
%    channel is an integer describing the channel number
%    metaData is a char array with whatever metaData you wish to include.
%       If it is over 64 characters, it will be truncated.
%    data is a 1xN numerical array, each element representing a single
%       ephys sample.
%    overwrite is an optional boolean flag indicating whether or not to
%       overwrite an existing file. Default false.
%
% This function is designed to write the provided 1D timeseries data and 
%   metadata to a binary file using the netCDF format. It is written with
%   a single channel of Intan ephys data in mind, and is designed to
%   produce files that can be read into electro_gui using the electro_gui
%   loader script egl_Intan_Bin.
%
% See also: readIntanNcFile, electro_gui, egl_Intan_Bin, 
%   convertIntanTxtToNc
%
% Version: 1.0
% Author:  Brian Kardon
% Email:   bmk27=cornell*org, brian*kardon=google*com
% Real_email = regexprep(Email,{'=','*'},{'@','.'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('overwrite', 'var')
    % Set default value for overwrite flag
    overwrite = false;
end

% Maximum length of metaData character string
metaDataLength = 64;

if length(metaData) < 64
    % Pad metaData with spaces to make it 64 char long
    metaData = [metaData, repmat(' ', [1, 64 - length(metaData)])];
elseif length(metaData) > 64
    % Warn user and truncate metaData to make it 64 char long
    warning('Warning, metaData has a max length of 64, and will be truncated');
    metaData = metaData(1:64);
end

if length(timeVector) == 1
    % This is a serial datenum. convert to 6-long time vector:
    timeVector = datevec(timeVector);
end

% Convert a 6-long time vector to 7-long by converting the factional
%   seconds into a separate microseconds entry.
if length(timeVector) == 6
    timeVector(7) = 1000000*timeVector(6)-1000000*floor(timeVector(6));
end


[path, filename, ext] = fileparts(filepath);
if isempty(ext)
    % No file extension found - add '.nc' on
    filepath = fullfile(path, [filename, '.nc']);
end

% Create blank netCDF file
if overwrite
    ncid = netcdf.create(filepath, 'CLOBBER');
else
    ncid = netcdf.create(filepath, 'NOCLOBBER');
end

try
    % Define variables that will go into file
    timeVectorDimID = netcdf.defDim(ncid,'p',7);
    timeVectorVarID = netcdf.defVar(ncid,'time','NC_INT',timeVectorDimID);
    deltaTVarID = netcdf.defVar(ncid,'dt','NC_DOUBLE',[]);
    channelVarID = netcdf.defVar(ncid,'chan','NC_INT',[]);
    metaDataDimID = netcdf.defDim(ncid,'c',metaDataLength);
    metaDataVarID = netcdf.defVar(ncid,'metaData','NC_CHAR',metaDataDimID);
    unlimDimId = netcdf.defDim(ncid,'t', netcdf.getConstant('NC_UNLIMITED'));
    dataVarID = netcdf.defVar(ncid,'data','NC_FLOAT',unlimDimId);
    netcdf.endDef(ncid);

    % Add data for each variable
    netcdf.putVar(ncid,timeVectorVarID,timeVector);
    netcdf.putVar(ncid,deltaTVarID,deltaT);
    netcdf.putVar(ncid,channelVarID,channel);
    netcdf.putVar(ncid,metaDataVarID,metaData);
    netcdf.putVar(ncid,dataVarID,0, length(data), data);
catch ME
    disp(getReport(ME));
end

% Close netCDF file.
netcdf.close(ncid)
