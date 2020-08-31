function [DataSetTransformed] = BiEx(DataSet,LookUpTable,PlotRange)
% Transforms FACS data (i.e. FACS_out{c,d}(:,fluorescentColor) into
% bi-exponential data according to: 
% 'A New "Logicle" Display Method Avoids Deceptive Effects of Logarithmic
% Scaling for Low Signals and Compensated Data. Parks et al. Cytometry Part
% A 69A:541-551 (2006)'
%
% Last modified: 05.05.2016
DataSetTransformed = interp1(LookUpTable,PlotRange,DataSet);
end