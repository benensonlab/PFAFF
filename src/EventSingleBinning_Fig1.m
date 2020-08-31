function binned_events = EventSingleBinning_Fig1(array,Bin_number)
% event-single-binning
%
% Bins data from .fcs files exported as .csv (through FlowJo) including ALL
% channels checked. Usually the columns are ordered as:
%
% xxx, xxx, BFP, Cerluean, Citrine, Cherry, xxx
%
% so only the columns 3-6 are needed. First, the data will be sorted along
% the transfection control to obtain the copy number bins. Only a range of
% intensities will be selected (2.5-97.5%). The bin sizes are determined
% according to the number of events within the intensity range. Afterwards
% the data is binned in equal (events) bins.
% 
% Second, each copy number bin will be binned again (events) by the inducer
% concentration (Cherry/PIT). The output of the function is a matrix of
% means from the binned data, number of BFP bins for each sample, number of
% Cherry bins for each BFP bin and the raw data divided into an array
% according to the bins of the format: [BFP,Cerulean,Citrine,Cherry]
%
% Last modified: 10.05.2016, CS

% Sort descending in the transfection control column
s_array = sortrows(array,-1); 

% Set intensity range
lb = prctile(s_array(:,1),0.1);
ub = prctile(s_array(:,1),99.9);

% Find events that are within the intensity range
Idx = find(s_array(:,1) >= lb & s_array(:,1) <= ub);

% Set the number of bins and compute the corresponding bin size
Bin_width = floor(length(Idx)/double(Bin_number));
bin_vector = 0:Bin_number;

% Select events, that fit all conditions
s_array = s_array(Idx,:);

% Pre-allocate
binned_events = cell(Bin_number,1);

for i = 1:Bin_number
    binned_events{i} = s_array(1+bin_vector(i)*Bin_width:bin_vector(i+1)*Bin_width,:);
end