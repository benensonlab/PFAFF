function binned_events = ValueSingleBinning_Fig1(array,Bin_number)
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
Idx = s_array(:,1) >= lb & s_array(:,1) <= ub;

% Select events, that fit all conditions
s_array = s_array(Idx,:);

% Define logscaled bin edges
BinEdges = linspace(lb,ub,Bin_number+1);

% Determine the bins for each entry
[~,~,Idx_Bin] = histcounts(s_array(:,1),BinEdges);

% Pre-allocate
binned_events = cell(Bin_number,1);

% Bin data
for i = 1:Bin_number
    binned_events{i} = s_array(i == Idx_Bin,:);
end