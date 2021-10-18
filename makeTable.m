clear
close all

% 
DataFolderName = '../../data/exampledata';                     % Edit to point to file folder
FileList = dir([DataFolderName, '/*/*.par']);        % list original datasets
nfile = length(FileList);
tableData = cell(nfile, 2);                         % use var name "tableData"

for subii = 1:nfile    
    subdir = FileList(subii).folder;
    subfile = FileList(subii).name;
    disp(['for ', subfile])

    % Add to data table
    tableData{subii, 1} = fullfile(subdir, subfile);        
    tableData{subii, 2} = 'Ready';
    
end

save('DataTable_template.mat', 'tableData');
disp('Done.')
