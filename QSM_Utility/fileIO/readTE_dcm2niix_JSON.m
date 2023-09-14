%% TE = readTE_dcm2niix_JSON(filenames)
%
% Input
% --------------
% filenames     : json file fullname(s)
%
% Output
% --------------
% TE            : echo times
%
% Description:
%
% Kwok-shing Chan @ DCCN
% k.chan@donders.ru.nl
% Date created: 13 June 2018
% Date last modified:
%
% modifed by X.L., using jsondecode

function [TE] = readTE_dcm2niix_JSON(filenames)

if ~iscell(filenames)
    filenames = {filenames};
end

TE =[];
for kf = 1:length(filenames)
    filename = filenames{kf};
    
    % using jsondecode
    text_json = fileread(filename);
    value_json = jsondecode(text_json);

    TE = [TE, value_json.EchoTime];     % in sec

end

TE = sort(TE,'ascend');

end
