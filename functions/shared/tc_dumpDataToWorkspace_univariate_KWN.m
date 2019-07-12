function tc_dumpDataToWorkspace_univariate_KWN
% Dumps data in the current fig (the current figure's 
% UserData) to the workspace. The variable name is 
% 'fig#Data', where # is the number of the current figure.
% The data is a struct.
%
% 06/22/04 ras.
tc = get(gcf,'UserData');
varname = sprintf('uvoxel_data');
assignin('base',varname,tc);
fprintf('Assigned plot data to struct %s in workspace.\n',varname);
return
