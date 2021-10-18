function R = Rmatrix_ang(rl, ap, fh, inv_flag, order)
% R = Rmatrix_ang(rl, ap, fh, inv_flag)

%% Header
% Purpose:
% make the rotation matrix according to the angulation parameter
% Assuming Tang = Tang = Trl*Tap*Tfh; 

% if inv_flag == 1
%     making the inverse rotation matrix i.e. Tang_inv =
%     Tfh_inv*Tap_inv*Trl_inv
% 
% INPUT
% rl, ap, fh are angulation parameters in degrees
% inv_flag tells makeing the rotation matrix or the inverse rotation matrix
% 
% OUTPUT 
% R is the 3X3 rotation matrix
%
% updated: X.L. 2013-05-21 added philips convention
% updated: X.L. 2016-12-07, Confirmed the definition with Philips
%               documentation, cross reference to DTI_gradient_table_creator.m in CATNAP

if nargin < 4
    inv_flag = 0;
    order = 1;
elseif nargin < 5
    order = 1;
end


T_AP = [ cosd(ap), 0, sind(ap); ...
         0, 1, 0; ...
        -sind(ap), 0, cosd(ap)];      % rotation about y, counter clockwise

T_FH = [cosd(fh), -sind(fh), 0;...
        sind(fh), cosd(fh), 0;...
        0, 0, 1];                   % rotation about z, counter clockwise

T_RL = [ 1, 0, 0; ...
        0, cosd(rl), -sind(rl); ...
        0, sind(rl), cosd(rl)];       % rotation about x, counter clockwise

if order == 1

    if inv_flag == 0    
        R = T_RL*T_AP*T_FH;             % Philips 1 % Correct, See also DTI_gradient_table_creator.m in CATNAP
    else
        R = T_FH'*T_AP'*T_RL';           
    end    

elseif order == 2
    disp('There may be some error, check TAng definition ...')
end