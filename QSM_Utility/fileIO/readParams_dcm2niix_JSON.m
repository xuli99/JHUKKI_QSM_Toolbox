%% Params = readParams_dcm2niix_JSON(filename, Params)
%
% Input
% --------------
% filenames     : json file fullname, Params
%
% Output
% --------------
% Params        : parameters contains 
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 June 2018
% Date last modified:
%
% modifed by X.L., using jsondecode

function Params = readParams_dcm2niix_JSON(filename, Params)

    text_json = fileread(filename);
    value_json = jsondecode(text_json);
    
    if isfield(value_json, 'MagneticFieldStrength')     % not in par/rec files, set in ParamsSetting file
        Params.B0 = value_json.MagneticFieldStrength;   % in Tesla
    end
    
    if isfield(value_json, 'RepetitionTime')
        Params.TR = value_json.RepetitionTime;          % in sec
    end
end
