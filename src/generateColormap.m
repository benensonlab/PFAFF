function [ColorMap,Colors] = generateColormap(varargin)
% Check for inputs, if there are any, scale the colormap accordingly
if nargin == 1
    Scale = varargin{1};
else
    Scale = 255;
end

% Colors are from Google's Material Design Palette
Red500 = [244,67,54];
lightOrange500 = [240,89,43];
Yellow500 = [254,234,57];
lightGreen500 = [139,194,72];
Green500 = [76,175,80];
Cyan500 = [9,188,211];
lightBlue500 = [56,164,220];
Blue500 = [33,150,243];
Indigo500 = [63,81,181];
deepPurple500 = [101,73,157];
White = [255,255,255];
Red900 = [183,28,28];

% Place all colors into one matrix
GradientColors = [Red500;
    %lightOrange500;
    Yellow500;
    lightGreen500;
    %Cyan500;
    lightBlue500;
    %Indigo500;
    %deepPurple500;
    White]./255;

% Place the colors along a vector (for simplicity, I use regular spacing,
% i.e. linspace)
x = [0,linspace(50,220,3),255];

% Generate actual colormap with interpolation between each consecutive
% color
ColorMap = interp1(x/255,GradientColors,linspace(1,0,Scale));
Colors = [Indigo500;Red500;Blue500;Green500;Red900;deepPurple500;Yellow500;Cyan500;lightBlue500;lightGreen500;lightOrange500]./255;

% % Testing
% I = linspace(0,1,255);
% imagesc(I(ones(1,10),:)');
% colormap(ColorMap)
