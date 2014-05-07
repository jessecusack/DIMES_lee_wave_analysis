function [varargout] = load_floats(floats)
% Compile the data from profiles into the classes Profile and EMApexFloat.
% Usage:
%   [EM4976a, EM4977a] = load_floats()
%
% NOTE: This function performs data quality control by checking that each
% profile has more than 100 pressure data points.
%
% There is also the option of interpolating GPS positions from a separate
% GPS file rather than using the profile GPS measure which is occasionally
% wrong.

if nargin < 1
    floats = {'4976a', '4977a'};
end

i = 1;
for float = floats
    
    if isunix
        root = ['/noc/users/jc3e13/data/EM-APEX/ema-' float{1} '-vel/'];
    else
        root = ['../../data/EM-APEX/ema-' float{1} '-vel/'];
    end % if
    
    dataFiles = dir([root 'ema-' float{1} '-*-vel.mat']);
    
    % In case the files were not added in numerical order.
%     datafiles = nestedSortStruct(datafiles, 'name', 1);
    
    profiles = [];
    
    % Transposition is necessary otherwise loop won't work.
    for dataFile = dataFiles'
        
        dataFilePath = fullfile(root, dataFile.name);
        
        data = load(dataFilePath);
        
        % Don't load floads with not enough data!!
        if length(data.Pctd) < 100; continue; end
        
        profiles = [profiles; Profile(data)];
        
    end % for over profiles
    
    varargout{i} = EMApexFloat(float{1}, profiles);
    i = i + 1;
    
end % for over floats

end % function