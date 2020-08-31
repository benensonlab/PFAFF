function [Mode,GoF,Condition] = GaussianFit_Secondary2Peaks(biex_data,N,BinCenter)
% Fits gaussians (gauss 3, gauss2 or gauss1) to biexponential
% distributions. Input data must look like:
% biex_data{c,d}{b,peaks},color(1)
%
% Last Update: 23.06.2016

% Initiate fitting options. Parameters MUST be positive
opt1 = fitoptions('gauss1','Lower',[0 0 0],'Upper',[Inf 4.5 4.5]); %,'StartPoint',[max(N),2.1,1]);
opt2 = fitoptions('gauss2','Lower',[0 0 0 0 0 0],'Upper',[Inf 4.5 4.5 Inf 4.5 4.5]); %,'StartPoint',[max(N),0.6,1,max(N),2.1,1]);
opt3 = fitoptions('gauss3','Lower',[0 0 0 0 0 0 0 0 0],'Upper',[Inf 4.5 4.5 Inf 4.5 4.5 Inf 4.5 4.5]); %,'StartPoint',[max(N),0.6,1,max(N),2.1,1,max(N),3.1,1]);

% Initiate empty/NaN vector
Empty = [NaN;NaN;NaN];

% Start with the fitting procedure

Gauss=0; % Records, which gauss fit is used
cond = 'NaN'; % Records, which fitting cascade condition is used

% Check if array is empty and if so, fill Mode with NaN and skip this loop
if isempty(biex_data)
    Mode = [Empty,Empty];
    GoF = NaN;
    Condition = NaN;
        
    % Compute the area under the curve
    Mode(4,1) = NaN;
    Mode(4,2) = NaN;
    return % go to the next loop iteration
end

% Start fitting cascade ---------------------------------------------------
% Fit gauss1 to data
[f1,gof1] = fit(BinCenter,N,'gauss1', opt1);

% CONDITION 1: Determine goodness of fit for gauss1, if worse than 97.5%,
% make gauss2
if gof1.rsquare > 0.95
    Mode = [[f1.a1;f1.b1;f1.c1],Empty];
    GoF = gof1.rsquare;
    
    Gauss=1;
    Condition = '1';
    
    % Compute the area under the curve
    Mode(4,1) = integral(@(x) Gauss_Curve(Mode(:,1),x),0,4.5);
    Mode(4,2) = NaN;
    return
else
    % Fit gauss2 to data
    [f2,gof2] = fit(BinCenter,N,'gauss2', opt2);
    Mode1 = [f2.a1;f2.b1;f2.c1]; % Peak b1
    Mode2 = [f2.a2;f2.b2;f2.c2]; % Peak b2
    gof = gof2.rsquare;
    
    Gauss=2;
    cond = '2a';
end

% Distance of peaks is smaller than 0.5: |b1-b2| < 0.5
if (Gauss>1) && (abs(Mode1(2)-Mode2(2)) < 0.5)
    % Fit gauss1 to data
    Mode = [[f1.a1;f1.b1;f1.c1],Empty];
    GoF = gof1.rsquare;
    
    Condition = '1b';
    Gauss=1;
    
    % Compute the area under the curve
    Mode(4,1) = integral(@(x) Gauss_Curve(Mode(:,1),x),0,4.5);
    Mode(4,2) = NaN;
    return
end

% Sort the highest peak into the first(!) column
if Gauss > 1
    if Mode1(2) >= Mode2(2)
        Mode = [Mode1,Mode2];
        GoF = gof;
    else
        Mode = [Mode2,Mode1];
        GoF = gof;
    end
end
Condition = cond;

% Compute the area under the curve
Mode(4,1) = integral(@(x) Gauss_Curve(Mode(:,1),x),0,4.5);
Mode(4,2) = integral(@(x) Gauss_Curve(Mode(:,2),x),0,4.5);

function [y] = Gauss_Curve(Pars,x)
% Compute the Gaussian for one peak.
a1=Pars(1,1);
b1=Pars(2,1);
c1=Pars(3,1);

y = a1*exp(-((x-b1)/c1).^2);