function [PhaseIn] = phase_unwrap_template(PhaseIn,TEs,TemplateEcho)
%[PhaseIn] = phase_unwrap_template(PhaseIn,TemplateEcho)
%   unwrap phase with template unwrapping
%   Reference: ROMEO paper, Dymerska et al, 2020, MRM
% 
% Xu Li
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

nEchoes = length(TEs);
UnwrapFlag = zeros(nEchoes, 1);  % 1: unwrapped

if TemplateEcho < 1 || TemplateEcho > nEchoes
    disp('Invalid selection of template echo... skip template unwrapping...')
    return;
end

PhaseIn(:,:,:,TemplateEcho) = phase_unwrap_path_mex(PhaseIn(:,:,:,TemplateEcho));
UnwrapFlag(TemplateEcho) = 1;

for echo_ind = [TemplateEcho:-1:1, TemplateEcho:nEchoes]
    if UnwrapFlag(echo_ind) < 1
        % still need to unwrap
        if echo_ind < TemplateEcho
            % unwrap echo_ind from echo_ind+1
            phaseTemplate = PhaseIn(:,:,:,echo_ind+1); TE_Template = TEs(echo_ind+1);
        elseif echo_ind > TemplateEcho    
            % unwrap echo_ind from echo_ind-1
            phaseTemplate = PhaseIn(:,:,:,echo_ind-1); TE_Template = TEs(echo_ind-1);
        end
        phaseToUnwrap = PhaseIn(:,:,:,echo_ind); TE_ToUnwrap = TEs(echo_ind);
        
        PhaseIn(:,:,:,echo_ind) = ...
            phaseToUnwrap - (2*pi)*(round((phaseToUnwrap - phaseTemplate*TE_ToUnwrap/TE_Template)/(2*pi)));

        UnwrapFlag(echo_ind) = 1; % mark
    end

end

end