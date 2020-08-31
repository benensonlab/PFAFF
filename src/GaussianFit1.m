function [Mode,GoF] = GaussianFit_Input(biex_data,N)
% Fits gaussians (gauss 3, gauss2 or gauss1) to biexponential
% distributions. Input data must look like: 
% Last Update: 23.06.2016
% Update: 03.09.2019, CS > differenct upper bounds for log scale

% Use global variables
global BinCenter

% Start with the fitting procedure
Gauss=0;
empty = [NaN;NaN;NaN];

% Check if array is empty and if so, skip. Fill Mode with NaN and skip this
% loop
if isempty(biex_data)
    Mode = [empty,empty];
    GoF = NaN;
    
    % Compute the area under the curve
    Mode(4,1) = NaN;
    Mode(4,2) = NaN;
    return
end;

% Start fitting cascade ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Fit gauss1 to data
options = fitoptions('gauss1','Lower',[0 0 0]); % Set constraints for parameters to be positive ,'Upper',[Inf 4.5 4.5],'StartPoint',[max(N),2.1,1]
[f1,gof1] = fit(BinCenter,N,'gauss1', options);

% Determine goodness of fit for gauss1, if worse than 97.5%, make gauss2
if gof1.rsquare > 0.95
    Mode = [[f1.a1;f1.b1;f1.c1],empty];
    GoF = gof1.rsquare;
    Gauss=1;
    
    % Compute the area under the curve
    Mode(4,1) = integral(@(x) Gauss_Curve(Mode(:,1),x),0,7);
    Mode(4,2) = NaN;
    return
else
    % Fit gauss2 to data
    options = fitoptions('gauss2','Lower',[0 0 0 0 0 0]); % Set constraints for parameters to be positive ,'Upper',[Inf 4.5 4.5 Inf 4.5 4.5],'StartPoint',[max(N),0.6,1,max(N),2.1,1]
    [f2,gof2] = fit(BinCenter,N,'gauss2', options);
    Mode1 = [f2.a1;f2.b1;f2.c1]; % Peak b1
    Mode2 = [f2.a2;f2.b2;f2.c2]; % Peak b2
    gof = gof2.rsquare;
    Gauss=2;
end;
% Determine goodness of fit for gauss 2; if worse than 95%, make gauss3
if gof2.rsquare > 0.99
    Mode1 = [f2.a1;f2.b1;f2.c1]; % Peak b1
    Mode2 = [f2.a2;f2.b2;f2.c2]; % Peak b2
    gof = gof2.rsquare;
    Gauss=2;
else
    % Fit gauss3 to data
    options = fitoptions('gauss3','Lower',[0 0 0 0 0 0 0 0 0]); % Set constraints for parameters to be positive ,'Upper',[Inf 4.5 4.5 Inf 4.5 4.5 Inf 4.5 4.5],'StartPoint',[max(N),0.6,1,max(N),2.1,1,max(N),3.1,1]
    [f3,gof3] = fit(BinCenter,N,'gauss3',options);
    
    % Sort gauss3 peaks by size (b)
    tmp_Mode1 = [f3.a1;f3.b1;f3.c1]; % Peak b1
    tmp_Mode2 = [f3.a2;f3.b2;f3.c2]; % Peak b2
    tmp_Mode3 = [f3.a3;f3.b3;f3.c3]; % Peak b3
    tmp = [tmp_Mode1,tmp_Mode2,tmp_Mode3];
    sort_tmp = sortrows(tmp',2)';
    
    % Decide which of the two peaks is higher (height(a)*SD(c)) and return
    % with that value as peak b1
    if sort_tmp(1,2)*sort_tmp(3,2)<sort_tmp(1,3)*sort_tmp(3,3)  %compare the "a" values of two highest "b"'s
        Mode1 = sort_tmp(:,3);
    else
        Mode1 = sort_tmp(:,2);
    end;
    Mode2 = sort_tmp(:,1); % 1 for lowest, 3 for highest
    gof = gof3.rsquare;
    Gauss=3;
    
    % Check, if the two highest peaks of gauss3 are too close and if, go
    % back to use the gauss2 fit
    if abs(sort_tmp(2,2)-sort_tmp(2,3))<0.3
        Mode1 = [f2.a1;f2.b1;f2.c1]; % Peak b1
        Mode2 = [f2.a2;f2.b2;f2.c2]; % Peak b2
        gof = gof2.rsquare;
        Gauss=2;
    end;
    
    % If the mean of the two tallest peaks from gauss3 is close to the tallest
    % peak of gauss2, use tall peak from gauss2 and smallest peak from
    % gauss3
    if (abs(mean([sort_tmp(2,2),sort_tmp(2,3)]) - max(f2.b1,f2.b2)) < 0.3)
        if f2.b1 > f2.b2 % tallest peak from gauss2
            Mode1 = [f2.a1;f2.b1;f2.c1];
        else
            Mode1 = [f2.a2;f2.b2;f2.c2];
        end;
        Mode2 = sort_tmp(:,1); % smallest peak from gauss3
        gof = gof2.rsquare;
        Gauss=2;
    end;
    
    clear tmp sort_tmp
end;

% Distance of peaks is smaller than 0.5: |b1-b2| < 0.5
if (Gauss==2) && (abs(Mode1(2)-Mode2(2)) < 0.75) %(Mode1(3)+Mode2(3))/2)
    % Fit gauss1 to data
    options = fitoptions('gauss1','Lower',[0 0 0]); % Set constraints for parameters to be positive ,'Upper',[Inf 4.5 4.5],'StartPoint',[max(N),2.1,1]
    [f1,gof1] = fit(BinCenter,N,'gauss1', options);
    Mode = [[f1.a1;f1.b1;f1.c1],empty];
    GoF = gof1.rsquare;
end;

% Sort the highest peak into the first(!) column
if Gauss > 1
    if Mode1(2) >= Mode2(2)
        Mode = [Mode1,Mode2];
        GoF = gof;
    else
        Mode = [Mode2,Mode1];
        GoF = gof;
    end;
end;

% Compute the area under the curve
Mode(4,1) = integral(@(x) Gauss_Curve(Mode(:,1),x),0,7);
Mode(4,2) = integral(@(x) Gauss_Curve(Mode(:,2),x),0,7);

function [y] = Gauss_Curve(Pars,x)
% Compute the Gaussian for one peak.
a1=Pars(1,1);
b1=Pars(2,1);
c1=Pars(3,1);

y = a1*exp(-((x-b1)/c1).^2);