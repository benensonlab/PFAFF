function OutData = RemoveSaturatedEvents(Data,ID,Saturation)
% Removes the saturated values from the data set. The saturation for our
% Fortessa is reached at 262143.
%
% (This script works only for case with up to 3 outputs! If there is a need
% for more, they have to be added. -- This issue was fixed.)
%
% 06.06.2017, Christoph Stelzer
% 27.03.2020, CS 

% switch ID.NrOut
%     case 1
%         Idx = Data(:,1) < Saturation & ...
%             Data(:,2) < Saturation & ...
%             Data(:,3) < Saturation;
%         
%     case 2
%         Idx = Data(:,1) < Saturation & ...
%             Data(:,2) < Saturation & ...
%             Data(:,3) < Saturation & ...
%             Data(:,4) < Saturation;
%         
%     case 3
%         Idx = Data(:,1) < Saturation & ...
%             Data(:,2) < Saturation & ...
%             Data(:,3) < Saturation & ...
%             Data(:,4) < Saturation & ...
%             Data(:,5) < Saturation;
% end

SzData = size(Data);
IDX = ones(SzData(1),1);

for Idx_Length = 1:SzData(2)
    Idx_notSaturated(:,Idx_Length) = Data(:,Idx_Length) < Saturation;    
    IDX = IDX & Idx_notSaturated(:,Idx_Length);
end
OutData = Data(IDX,:);