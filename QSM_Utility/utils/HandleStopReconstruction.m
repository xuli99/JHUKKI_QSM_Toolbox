%% Author: Jiri van Bergen
% Affiliation: Radiology @ JHU - Kirby Center
% Contact via xuli@mri.jhu.edu

% IF USER PRESSES THE X BUTTON
% ALL THE STOPPING STUFF 

if(hasCanceled) 
    % Close it all
    multiWaitbar('CloseAll');
    
    % Enable buttons
    set([handles.ButtonEditEchoes handles.VarFSLThres handles.VarB0 handles.VarRadiusDisk handles.VarMaskEchoes handles.VarFSL handles.VarBgRemoval handles.VarSHARPradius handles.VarQSMSolver], 'Enable', 'On');

    % Stop
    error('Reconstruction stopped by user'); 
end