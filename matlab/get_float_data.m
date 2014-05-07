function [floatStruct] = getFloatData(floatID)
% Extract properties for a single float from an allprof##.mat file. 

% Make this an argument of the function for added generality. Not necessary
% for the moment. 
if isunix
    AP = load('/noc/users/jc3e13/data/EM-APEX/allprofs11.mat');
else
    AP = load('../../data/EM-APEX/allprofs11.mat');
end
% Check if float ID entered correctly.
if ismember(floatID, AP.flid)
    inFloat = AP.flid == floatID;
else
    error('Specified float cannot be found.')
end

% Geographic limits applies to all profiles in file and will cause problems
% in loop. Remove.
AP = rmfield(AP, 'ar');

% Use some known fields to get float property sizes.
NQ = size(AP.P);
MQ = size(AP.U);
PQ = size(AP.P_ca);
Q = size(AP.flid);

floatStruct = struct();

for fieldname = fieldnames(AP)'
    % Get name out of cell.
    fn = fieldname{1};
    fieldsize = size(AP.(fn));
    % Different size array will need slightly different indexes.
    if fieldsize == Q
        floatStruct.(fn) = AP.(fn)(inFloat);
    elseif fieldsize == PQ | fieldsize == MQ | fieldsize == NQ
        floatStruct.(fn) = AP.(fn)(:, inFloat);
    end
end