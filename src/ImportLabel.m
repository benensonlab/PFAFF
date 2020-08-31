function [ID] = ImportLabel(FILE)
% Open Label File
LabelfileID = fopen(FILE);
% Read file and store the labels and column IDs
LabelID = textscan(LabelfileID,'%s %u','CommentStyle','%');
% Close 'Description' file
fclose(LabelfileID);

ID.Circuits = LabelID{2}(1);  % Number of Circuits
ID.InputLvl = LabelID{2}(2);  % Number of Input Levels
ID.TMBins = LabelID{2}(3);    % Number of TransfectionMarker Bins

ID.colTM = LabelID{2}(4); % Column TransfectionMarker
ID.colIn = LabelID{2}(5); % Column Input
ID.colO1 = LabelID{2}(6); % Column Output_1

% Determine the number of outputs
ID.NrOut = size(LabelID{1});

switch ID.NrOut(1)
    case 7
        ID.colO2 = LabelID{2}(7); % Column Output_2
        
    case 8
        ID.colO2 = LabelID{2}(7); % Column Output_2
        ID.colO3 = LabelID{2}(8); % Column Output_3
end;