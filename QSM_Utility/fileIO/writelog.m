function writelog(LogFile, LogText, perm)
% writing to log file for cluster version
if nargin < 3
    perm = 'a';
end
fID = fopen(LogFile, perm);
fprintf(fID, LogText);
fclose(fID);