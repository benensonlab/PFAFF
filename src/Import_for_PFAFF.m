function [DataStructure,ID] = Import_for_PFAFF
% This function imports .csv-files from a folder into a structure usable
% for PFAFF. The arguments define the:
%
% NrReplicates: Number of replicates 
% NrCircuits: Number of circuits 
% NrInputs: Number of inputs
%
% The output has the form:
% <DataStructure>.<Replicate>.<FlowData>{Idx_Circuit,Idx_Input}
%
% 19.04.2020, CS

% Fetch all files from the folder
files = dir('*.csv');
% For each file in the folder do
for Idx_files = 1:length(files)
    % Split the filename at '_' characters
    filename{Idx_files,:}= regexp(files(Idx_files).name,'_','split');
    % Remove .csv tag from the filenames
    filename{Idx_files,1}{end}(end-3:end) = [];
    
    % Extract all input modulation levels
    Input(Idx_files,:) = str2double(filename{Idx_files}{3}(regexp(filename{Idx_files}{3},'\d')));
    % Extract the input modulation level names
    Inputs{Idx_files,:} = filename{Idx_files}{3};
    % Extract the circuit names
    Circuits{Idx_files,:} = filename{Idx_files}{2};
    % Extract the replicate names
    Replicates{Idx_files,:} = filename{Idx_files}{4};
end

% Collect all unique names in the structure
ID.Circuit = unique(Circuits);
ID.Input = unique(Inputs);
ID.Replicate = unique(Replicates);

% Loop through all filenames and assign each filename to the correct
% position in DataStructure
for Idx_files = 1:length(files)
    for r = 1:length(ID.Replicate)
        for f = 1:length(ID.Circuit)
            for z = 1:length(ID.Input)
                if strcmp(filename{Idx_files}{2},ID.Circuit{f}) && ...
                        strcmp(filename{Idx_files}{3},ID.Input{z}) && ...
                        strcmp(filename{Idx_files}{4},ID.Replicate{r})
                    DataStructure.(sprintf('Rep%i',r)).FlowData{f,z} = csvread(files(Idx_files).name,1,0);
                end
            end
        end
    end
end