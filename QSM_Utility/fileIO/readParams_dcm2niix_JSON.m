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
% modifed by X.L.

function Params = readParams_dcm2niix_JSON(filename, Params)

    % read mode
    fid = fopen(filename,'r');
    % read file line by line
    line = fgetl(fid);
    
    % start reading lines
    while ischar(line)

        if contains(line, 'MagneticFieldStrength')
            Params.B0 = get_str(line);
        end

        if contains(line, 'RepetitionTime')
            Params.TR = get_str(line);
        end
        
        % start the next line
        line = fgetl(fid);
    end

    fclose(fid);

end

%% Get value of the tag 
function str=get_str(list_info)
    % find chars ': '
    k_b = strfind(list_info,': ');
    % account for two chars ':' and ' '
    str=str2double(list_info(k_b(1)+2:end));
end