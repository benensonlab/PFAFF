function [BackTransformed] = invBiEx(DataSetTransformed,varargin)
% Transforms biexponentially transformed data (BiEx.m) back into FACS
% units. Source Paper: 'A New "Logicle" Display Method Avoids Deceptive
% Effects of Logarithmic Scaling for Low Signals and Compensated Data.
% Parks et al. Cytometry Part A 69A:541-551 (2006)'
%
% Last modified: 15.11.2016

% Determine the number of input variables 
numvarargs = length(varargin);

% Show error message for more than 4 input variables
if numvarargs > 4
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 4 optional inputs');
end

% Set degault values for optinal inputs
optargs = {4.5 2 2^18 0.401};

% Optional change of default values
optargs(1:numvarargs) = varargin;
[M,p,T,W] = optargs{:};
% M = 4.5;    % length of display
% p = 2;      % parameter for compactness, but p & W are one adjustable parameter
% T = 2^18;   % max data value
% W = 0.401;  % Strength and range of linearization around 0

% Determine the size of the input array
Input_Size = size(DataSetTransformed);

% Index for values under case 1 (i.e. input is larger or equal than W)
Idx = DataSetTransformed >= W;
% Back transformation of case 1
delta = DataSetTransformed(Idx)-W;
BackTransformed(Idx) = T.*10^(-(M-W)).*(10.^(delta)-p.^2.*10.^(-(delta)./p)+p.^2-1);
% Back transformation of case 2 (= else)
delta = W-DataSetTransformed(~Idx);
BackTransformed(~Idx) = -T.*10.^(-(M-W)).*(10.^(delta)-p.^2.*10.^(-(delta)./p)+p.^2-1);
% Reshape the back transformed to the input size array
BackTransformed = reshape(BackTransformed,Input_Size);