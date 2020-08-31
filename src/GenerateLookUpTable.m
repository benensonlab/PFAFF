function [LookUpTable,PlotRange] = GenerateLookUpTable(varargin)
% Generate LookUpTable for a given M (length of display), p (parameter for
% compactness, but p & W are one adjustable parameter), T (max data value),
% W (Strength and range of linearization around 0).

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

% Compute the LookUpTable
PlotRange = 0:0.01:4.5;
j=1;
for x = PlotRange
    if x >= W
        delta = x-W;
        LookUpTable(j) = T*10^(-(M-W))*(10^(delta)-p^2*10^(-(delta)/p)+p^2-1);
    else
        delta = W-x;
        LookUpTable(j) = -T*10^(-(M-W))*(10^(delta)-p^2*10^(-(delta)/p)+p^2-1);
    end
    j = j+1;
end
end