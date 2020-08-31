function [DataSet,Data,ID] = PFAFF(Data,NumberOfBins,colTM,colIn,varargin)
%PFAFF maps the input to the circuits output of flow cytometry data.
%
%DESCRIPTION
% Genetic circuits at various input levels are analyzed via flow cytometry
% and the corresponding input/output relations are extracted from these
% data sets. This function wraps all necessary modules (workspaces,
% sub-functions, etc.) and saves the output in the folder './PFAFF_Results'
%
%SYNTAX
% [DataSet,DATA,ID] = PFAFF(Data,NumberOfBins,ColumnTM,ColumnInput,ColumnOutput1,ColumnOutput2,...)
% 
%INPUT
% Data is a structure that contains multiple replicates (e.g. DATA.Rep1).
%   Each replicate is itself a structure containing the flow data cell array,
%   in which each row represents a genetic circuit and each column an input
%   level (e.g. DATA.Rep1.Flow{CircuitA,Input1}). Replicates should have the
%   same number of circuits and inputs or at least should be placed in the
%   same array position.
% NumberOfBins is the number of bins for the transfection marker
% ColumnTM is the column position of the transfection marker in the data
% set
% ColumnInput is the column position of the input in the data set
% ColumnOutput1 is the column position of the output1 in the data set
% ColumnOutput2 is the column position of the output2 in the data set
% ...
%
%OUTPUT
% DataSet a structure containing all extracted mode positions for each
% circuit and bin
%   A{Idx_Color}{Idx_Circuit,Idx_Bin} are the mode position's peak heights
%   S{Idx_Color}{Idx_Circuit,Idx_Bin} are the output's modes with the
%   highest weights
%   W{Idx_Color}{Idx_Circuit,Idx_Bin} are the output's mode weights
%   X{Idx_Circuit,Idx_Bin} are the input's mode positions 
%   Y{Idx_Color}{Idx_Circuit,Idx_Bin} are the output's mode positions
%
% DATA a structure contains all arrays of intermediate computation steps
% ID a structure containing all label information
%
%EXAMPLE
% [Result,~,~] = PFAFF(bb,10,3,6,4,5);
%
% 05.06.2017, Christoph Stelzer

%% Start timer ********************************************************** %
tic;

% Add source folder to path
addpath('src');

% Bi-exponential transformation
biex.M = 4.5; % length of display
biex.p = 2; % parameter for compactness, but p & W are one adjustable parameter
biex.T = 2^18; % max data value
biex.W = 0.401; % Strength and range of linearization around 0
% Generate the LookUpTable and Plotrange for the given biex parameter set
[LookUpTable,PlotRange] = GenerateLookUpTable(biex.M,biex.p,biex.T,biex.W);
% Save the 
save('LookUpTable.mat','LookUpTable','PlotRange');

% Import labels and information about the input data into the structure ID
for Idx_arg = 1:length(varargin)
    % Look through all variable arguments for the Name-Value pair
    if strcmpi(varargin{Idx_arg},'Labels')
        % If the Name-Value pair is found, assign it
        ID = varargin{Idx_arg+1}; 
        % Delete the entries from the input arguments
        varargin{Idx_arg+1}=[]; 
        varargin{Idx_arg}=[];
        FlagID = true;
    end
end

ID.Replicates = length(fieldnames(Data)); % Number of replicates
ID.ReplicateField = fieldnames(Data); % Name of variables on replicate level
ID.DataVariable = fieldnames(Data.(ID.ReplicateField{1})); % Name of variables on the data level
ID.DataSize = size(Data.(ID.ReplicateField{1}).(ID.DataVariable{1})); % Number of rows and columns of input data
ID.Circuits = ID.DataSize(1); % Number of Circuits
ID.InputLvl = ID.DataSize(2); % Number of Input Levels
ID.TMBins = NumberOfBins; % Number of TransfectionMarker Bins
ID.colTM = colTM; % Column TransfectionMarker
ID.colIn = colIn; % Column Input
Idx_Out = 1;
if exist('FlagID','var')
    RequiredInArg = 6;
else
    RequiredInArg = 4;
    for Idx_Circuit = 1:ID.Circuits
        ID.Circuit{Idx_Circuit} = char('A'+(Idx_Circuit)-1);
    end
end
while (Idx_Out <= nargin - RequiredInArg)
    ID.colOut(Idx_Out) = varargin{Idx_Out}; % Column Output_x
    ID.Label{Idx_Out+2} = sprintf('Output_%i',Idx_Out); % Label Output_x
    Idx_Out = Idx_Out + 1;
end
ID.NrOut = nargin - RequiredInArg; % Number of outputs
ID.StructName = inputname(1); % Variable name of data
% Additional labels
% ID.DOX = [0,1,3.5,12,42,144,500,1500];
ID.Flow_scale = BiEx([0,100,1000,10000,100000],LookUpTable,PlotRange);
ID.Label{1} = 'TransfectionMarker';
ID.Label{2} = 'Input';
ID.Peak = {'high','low'};

% Create result path and folder
check_fldr = exist(sprintf('PFAFF_Result_%s',ID.StructName),'dir');
if check_fldr == 7
    t = datetime('now');
    folder = fullfile(pwd,sprintf('PFAFF_Result_%s_%s',ID.StructName,datestr(t,30)));
else
    folder = fullfile(pwd,sprintf('PFAFF_Result_%s',ID.StructName));
end
mkdir(folder);

% Generate colors and colormap
[MaterialMap,MaterialColors] = generateColormap;

% Preallocate variables
ub = cell(ID.Circuits,ID.InputLvl,ID.TMBins);
lb = cell(ID.Circuits,ID.InputLvl,ID.TMBins);
Idx_Slice = cell(ID.Circuits,ID.InputLvl,ID.TMBins);
Flow_pl = cell(ID.Circuits,ID.InputLvl);
RepSize = cell(ID.Replicates,1);

% If the data set contains replicates, check if they have the same size
if ID.Replicates > 1
    for Idx_Rep = 1:ID.Replicates
        RepSize{Idx_Rep} = size(Data.(ID.ReplicateField{Idx_Rep}).(ID.DataVariable{1}));
    end
    msg = 'Error. Replicate arrays have different sizes.';
    if ~isequal(RepSize{:})
        error(msg);
    end
end

% Set global variables
global BinEdge BinCenter
BinEdge = linspace(0,4.5,101);
HalfBin = (BinEdge(2) - BinEdge(1))/2;
BinCenter = BinEdge(1:100)' + HalfBin;

%% Binning data according to the transfection marker
for Idx_Rep = 1:ID.Replicates
    % Fetch data variable name
    DataVariable = fieldnames(Data.(ID.ReplicateField{Idx_Rep}));
    
    for Idx_Circuit = 1:ID.Circuits
        % Initialize tmp variables
        PeakValue = zeros(8,1);
        for Idx_Input = 1:ID.InputLvl
            % Reorder input array
            Data.(ID.ReplicateField{Idx_Rep}).DataArray{Idx_Circuit,Idx_Input} = horzcat(...
                Data.(ID.ReplicateField{Idx_Rep}).(DataVariable{1}){Idx_Circuit,Idx_Input}(:,ID.colTM),...
                Data.(ID.ReplicateField{Idx_Rep}).(DataVariable{1}){Idx_Circuit,Idx_Input}(:,ID.colIn),...
                Data.(ID.ReplicateField{Idx_Rep}).(DataVariable{1}){Idx_Circuit,Idx_Input}(:,ID.colOut));
            
            % Set saturation threshold
            Saturation = 2^18-1;
            
            % Remove saturated values
            Data.(ID.ReplicateField{Idx_Rep}).DataArray{Idx_Circuit,Idx_Input} = RemoveSaturatedEvents(Data.(ID.ReplicateField{Idx_Rep}).DataArray{Idx_Circuit,Idx_Input},ID,Saturation);
            
            % Bin data %original
            Data.(ID.ReplicateField{Idx_Rep}).binned_events{Idx_Circuit,Idx_Input} = EventSingleBinning(Data.(ID.ReplicateField{Idx_Rep}).DataArray{Idx_Circuit,Idx_Input},ID.TMBins);
            % Bin data %according to dPFAFF
            % [Data.(ID.ReplicateField{Idx_Rep}).binned_events{Idx_Circuit,Idx_Input},PeakValue(Idx_Input)] = PeakBinning(Data.(ID.ReplicateField{Idx_Rep}).DataArray{Idx_Circuit,Idx_Input},3,Idx_Input,PeakValue,LookUpTable,PlotRange);

            % Transform binned data biexponentially
            for Idx_Bin = 1:ID.TMBins
                Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin} = BiEx(Data.(ID.ReplicateField{Idx_Rep}).binned_events{Idx_Circuit,Idx_Input}{Idx_Bin},LookUpTable,PlotRange);
            end
            
            % Concatenate biex data of all input levels for each bin
            for Idx_Bin = 1:ID.TMBins
                Data.(ID.ReplicateField{Idx_Rep}).comp_biex{Idx_Circuit,Idx_Bin} = vertcat(Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,:,Idx_Bin});
            end
        end
%         ID.PeakValue(Idx_Circuit,:) = PeakValue;
        clear PeakValue
    end
    
    %% Fit Input
    % Create output folder(path)
    subfolder_Input = fullfile(folder,ID.ReplicateField{Idx_Rep},'InputPlots');
    mkdir(subfolder_Input);
    
    Idx_Color = 2; % Input column
    Peak = 1; % There is only one 'InputPeak' for the input
    
    for Idx_Circuit = 1:ID.Circuits
        for Idx_Input = 1:ID.InputLvl
            for Idx_Bin = 1:ID.TMBins
                % Store histcounts
                Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color) = histcounts(Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin}(:,Idx_Color),BinEdge);
                
                % Extract input modes
                Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin} = GaussianFit_Input(Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin}(:,Idx_Color),Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin}(:,Idx_Color));
                
                % Back-transform (from biexponential space) into linear space
                Data.(ID.ReplicateField{Idx_Rep}).invMode_Input_dox{Idx_Circuit,Idx_Bin}(Idx_Input,:) = invBiEx(Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,:));
            end
        end
    end
    
    %% Plot distributions and peaks to verify strategy
    for Idx_Bin = 1:ID.TMBins
        fig = figure;
        for Idx_Input = 1:ID.InputLvl
            for Idx_Circuit = 1:ID.Circuits
                % Extract mode's x, y positions
                clear a1 a2 b1 b2
                a1 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(1,1);
                b1 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,1);
                a2 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(1,2);
                b2 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,2);
                % Plot results
                subplot(double(ID.Circuits),double(ID.InputLvl),double(Idx_Input+(Idx_Circuit-1)*ID.InputLvl))
                plot(BinCenter,Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin,1}(:,Idx_Color),'.k');
                hold on;
                % Label the peak modes
                text(b1,a1,'+b1','color','r');
                text(b2,a2,'+b2','color','r');
                hold off
                % Annotate the plot
                title(sprintf('%s%i, InputLvl %i',ID.Circuit{Idx_Circuit},Idx_Bin,Idx_Input));
                xlabel(sprintf('%s',ID.Label{Idx_Color}));
                ylabel('Counts');
                xlim([0 4.5]);
                verticalLim = ylim;
                ylim([0 verticalLim(2)]);
            end
        end
        % Make plot screen size and save it as .png
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        fig.PaperPositionMode = 'auto';
%         fig.Visible = 'on';
%         % Save figure as .png-file
%         print(fullfile(subfolder_Input,sprintf('Bin%i_Mode_%s',Idx_Bin,ID.Label{Idx_Color})),'-dpng');
        % Save figure as .fig-file
        savefig(fullfile(subfolder_Input,sprintf('Bin%i_Mode_%s_Rep%i.fig',Idx_Bin,ID.Label{Idx_Color},Idx_Rep)));
        close all hidden;
    end
    
    %% Fit Output
    % Create directories for outputs
    Idx_Out = 1;
    while Idx_Out <= ID.NrOut
        subfolder_Output{Idx_Out} = fullfile(folder,ID.ReplicateField{Idx_Rep},sprintf('Output%i_Plots',Idx_Out));
        mkdir(subfolder_Output{Idx_Out});
        Idx_Out = Idx_Out + 1;
    end
    
    for Idx_Circuit = 1:ID.Circuits
        for Idx_Input = 1:ID.InputLvl
            for Idx_Bin = 1:ID.TMBins
                % Set boundaries around input peak for secondary binning
                lb{Idx_Circuit,Idx_Input,Idx_Bin} = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,:) - 0.1;
                ub{Idx_Circuit,Idx_Input,Idx_Bin} = Data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,:) + 0.1;
                
                % Extract all events within this boundaries
                for Peak = 1:2
                    % Cut slice around input peak
                    Idx_Slice{Idx_Circuit,Idx_Input,Idx_Bin} = Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin}(:,2) >= lb{Idx_Circuit,Idx_Input,Idx_Bin}(Peak) & Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin}(:,2) <= ub{Idx_Circuit,Idx_Input,Idx_Bin}(Peak);
                    
                    % Extract all events within in this slice
                    Data.(ID.ReplicateField{Idx_Rep}).biex_data_Chypeak{Idx_Circuit,Idx_Input,Idx_Bin,Peak} = Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,Idx_Bin}(Idx_Slice{Idx_Circuit,Idx_Input,Idx_Bin},:);
                    
                    % Building histogram data
                    for Idx_Color = 3:ID.NrOut+2
                        % Compute histcounts
                        Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color) = histcounts(Data.(ID.ReplicateField{Idx_Rep}).biex_data_Chypeak{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color),BinEdge);
                        
                        % Preallocation
                        tmpSpline = cell(Idx_Circuit,Idx_Input,Idx_Bin,2,5);
                        SplineXY = cell(Idx_Circuit,Idx_Input,Idx_Bin,2,5);
                        
                        % Smooth data by using Splines
                        tmpSpline{Idx_Circuit,Idx_Input,Idx_Bin,Peak,Idx_Color} = csaps(BinCenter,Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color),0.999);
                        SplineXY{Idx_Circuit,Idx_Input,Idx_Bin,Peak,Idx_Color} = fnplt(tmpSpline{Idx_Circuit,Idx_Input,Idx_Bin,Peak,Idx_Color});
                        
                        Data.(ID.ReplicateField{Idx_Rep}).SplineX{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color) = SplineXY{Idx_Circuit,Idx_Input,Idx_Bin,Peak,Idx_Color}(1,:);
                        Data.(ID.ReplicateField{Idx_Rep}).SplineY{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color) = SplineXY{Idx_Circuit,Idx_Input,Idx_Bin,Peak,Idx_Color}(2,:);
                        clear tmpSpline
                    end
                    Idx_Out = 1;
                    while Idx_Out <= ID.NrOut
                        % Gaussian Fit of Output
                        Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{Peak,Idx_Circuit,Idx_Input,Idx_Bin} = GaussianFit_Secondary(...
                            Data.(ID.ReplicateField{Idx_Rep}).biex_data_Chypeak{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Out+2),...
                            Data.(ID.ReplicateField{Idx_Rep}).SplineY{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Out+2),...
                            Data.(ID.ReplicateField{Idx_Rep}).SplineX{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Out+2));
                        
                        % Back-transform (from biexponential space) into linear space
                        Data.(ID.ReplicateField{Idx_Rep}).invMode_Output_dox{Idx_Out}{Peak,Idx_Circuit,Idx_Bin}(Idx_Input,:) = invBiEx(Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{Peak,Idx_Circuit,Idx_Input,Idx_Bin}(2,:));
                        Idx_Out = Idx_Out + 1;
                    end
                end
            end
        end
    end
    
    %% Generate histogram of biex data, plot and indicate the peaks from fitting
    for Idx_Color = 3:ID.NrOut+2
        for Peak = 1:2
            for Idx_Bin = 1:ID.TMBins
                fig = figure;
                for Idx_Input = 1:ID.InputLvl
                    for Idx_Circuit = 1:ID.Circuits
                        % Extract mode's x, y positions
                        clear a1 a2 b1 b2
                        a1 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Color-2}{Peak,Idx_Circuit,Idx_Input,Idx_Bin}(1,1);
                        b1 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Color-2}{Peak,Idx_Circuit,Idx_Input,Idx_Bin}(2,1);
                        a2 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Color-2}{Peak,Idx_Circuit,Idx_Input,Idx_Bin}(1,2);
                        b2 = Data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Color-2}{Peak,Idx_Circuit,Idx_Input,Idx_Bin}(2,2);
                        
                        % Plot results
                        subplot(double(ID.Circuits),double(ID.InputLvl),double(Idx_Input+(Idx_Circuit-1)*ID.InputLvl))
                        plot(BinCenter,Data.(ID.ReplicateField{Idx_Rep}).N{Idx_Circuit,Idx_Input,Idx_Bin,Peak}(:,Idx_Color),'.k');
                        hold on;
                        % Label the peak modes
                        text(b1,a1,'+b1','color','r');
                        text(b2,a2,'+b2','color','r');
                        hold off
                        % Annotate the plot
                        title(sprintf('%s%i, InputLvl %i',ID.Circuit{Idx_Circuit},Idx_Bin,Idx_Input));
                        xlabel(sprintf('%s',ID.Label{Idx_Color}));
                        ylabel('Counts');
                        xlim([0 4.5]);
                        verticalLim = ylim;
                        ylim([0 verticalLim(2)]);
                    end
                end
                % Make plot screen size and save it as .png
                fig.Units = 'normalized';
                fig.OuterPosition = [0 0 1 1];
                fig.PaperPositionMode = 'auto';
%                 fig.Visible = 'on';
%                 % Save figure as .png-file
%                 print(fullfile(subfolder_Output{Idx_Color-2},sprintf('Bin%i_Mode_Output%i_InputPeak%i',Idx_Bin,Idx_Color-2,Peak)),'-dpng');
                % Save figure as .fig-file
                savefig(fullfile(subfolder_Output{Idx_Color-2},sprintf('Bin%i_Mode_Output%i_Input_mode_%s_Rep%i.fig',Idx_Bin,Idx_Color-2,ID.Peak{Peak},Idx_Rep)));
                close all hidden;
            end
        end
    end
    
    %% Flow plots
    % Output folder
    folder_flow = fullfile(folder,ID.ReplicateField{Idx_Rep},'Flow_Plots');
    mkdir(folder_flow);
    
    % Select colors
    for XColor = 2:ID.NrOut+2
        for YColor = 3:ID.NrOut+2
            % Omit plots that show the same color on x and y-axis
            if (XColor == YColor) || (XColor > YColor)
                continue
            end
            
            % Create FACS Plots
            fig = figure;
            for Idx_Circuit = 1:ID.Circuits
                for Idx_Input = 1:ID.InputLvl
                    % Concatenate bins
                    Flow_pl{Idx_Circuit,Idx_Input} = vertcat(Data.(ID.ReplicateField{Idx_Rep}).biex_data{Idx_Circuit,Idx_Input,:});
                    
                    % Plot
                    subplot(double(ID.Circuits),double(ID.InputLvl),double(Idx_Input+(Idx_Circuit-1)*ID.InputLvl))
                    cloudPlot(Flow_pl{Idx_Circuit,Idx_Input}(:,XColor),Flow_pl{Idx_Circuit,Idx_Input}(:,YColor),[0 4.5 0 4.5]);
                    xlabel(sprintf('%s',ID.Label{XColor}));
                    ylabel(sprintf('%s',ID.Label{YColor}));
                    title(sprintf('%s, Inputlvl %i',ID.Circuit{Idx_Circuit},Idx_Input));
                    grid on;
                    ax = gca;
                    ax.XTick = ID.Flow_scale; ax.XTickLabel = {'0','10^2','10^3','10^4','10^5'};
                    ax.YTick = ID.Flow_scale; ax.YTickLabel = {'0','10^2','10^3','10^4','10^5'};
                    
                    % Change to pretty colors (Material Design Colorpalette)
                    colormap(MaterialMap);
                    caxis([0 2000]);
                end
            end
            
            % Make plot screen size and save it as .png
            fig.Units = 'normalized';
            fig.OuterPosition = [0 0 1.08/1.92 1];
            fig.PaperPositionMode = 'auto';
            fig.Visible = 'on';
%           % Save figure as .png-file
            print(fullfile(folder_flow,sprintf('FlowPlot_%s_v_%s.png',ID.Label{XColor},ID.Label{YColor})),'-dpng');
            % Save figure as .fig-file
            savefig(fullfile(folder_flow,sprintf('FlowPlot_%s_v_%s_Rep%i.fig',ID.Label{XColor},ID.Label{YColor},Idx_Rep)));
            close all hidden
        end
    end
    
    % Clean up struture
    Data.(ID.ReplicateField{Idx_Rep}) = rmfield(Data.(ID.ReplicateField{Idx_Rep}),{'SplineX','SplineY','invMode_Input_dox','invMode_Output_dox'});
end

%% Re-structure the results
DataSet = Modes2Array(Data,ID);

% If there are replicates, compute the means and SD
if ID.Replicates > 1
    DataSet = MakeMeans(DataSet,ID);
end

%% Bubbleplots of modes

% Output folder
folder_bubble = fullfile(folder,'IO_Plots');
mkdir(folder_bubble);

for Idx_Color = 3:ID.NrOut+2
    for Idx_Bin = 1:ID.TMBins
        fig = figure;
        for Idx_Circuit = 1:ID.Circuits
            % Plot flow recordings
            subplot(double(ceil(sqrt(ID.Circuits))),double(ceil(sqrt(ID.Circuits))),double(Idx_Circuit))
            scatter(DataSet.X{Idx_Circuit,Idx_Bin}(:),DataSet.Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:),100*DataSet.W{Idx_Color}{Idx_Circuit,Idx_Bin}(:),'markerfacecolor',MaterialColors(Idx_Circuit,:),'markeredgecolor','k');
            % Annotate plot
            title(sprintf('%s Bin%i',ID.Circuit{Idx_Circuit},Idx_Bin));
            xlabel(ID.Label{2});
            ylabel(ID.Label{Idx_Color});
            ax = gca;
            ax.XScale = 'log'; xlim([1 2e5]);
            ax.YScale = 'log'; ylim([1 2e5]);
            grid on
        end
        
        % Make plot screen size and save it as .png
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        fig.PaperPositionMode = 'auto';
        fig.Visible = 'on';
%       % Save figure as .png-file
        print(fullfile(folder_bubble,sprintf('Input_%s_Bin%i.png',ID.Label{Idx_Color},Idx_Bin)),'-dpng');
        % Save figure as .fig-file
        savefig(fullfile(folder_bubble,sprintf('Input_%s_Bin%i.fig',ID.Label{Idx_Color},Idx_Bin)));
        close all hidden;
    end
end

% If there are replicates, plot the mean modes and with their SDs
if ID.Replicates > 1
    % Output folder
    folder_bubble = fullfile(folder,'IO_Mean');
    mkdir(folder_bubble);
    
    for Idx_Color = 3:ID.NrOut+2
        for Idx_Bin = 1:ID.TMBins
            fig = figure;
            for Idx_Circuit = 1:ID.Circuits
                % Plot flow recordings
                subplot(double(ceil(sqrt(ID.Circuits))),double(ceil(sqrt(ID.Circuits))),double(Idx_Circuit))
                scatter(DataSet.X_bar{Idx_Circuit,Idx_Bin}(:),DataSet.Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:),100*DataSet.W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:),'markerfacecolor',MaterialColors(Idx_Circuit,:),'markeredgecolor','k');
                hold on
                % Errorbars x-axis
                line([DataSet.X_SD_lb{Idx_Circuit,Idx_Bin}(:),DataSet.X_SD_ub{Idx_Circuit,Idx_Bin}(:)]',[DataSet.Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:),DataSet.Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:)]','color','k','linewidth',0.1);
                % Errorbars y-axis
                line([DataSet.X_bar{Idx_Circuit,Idx_Bin}(:),DataSet.X_bar{Idx_Circuit,Idx_Bin}(:)]',[DataSet.Y_SD_lb{Idx_Color}{Idx_Circuit,Idx_Bin}(:),DataSet.Y_SD_ub{Idx_Color}{Idx_Circuit,Idx_Bin}(:)]','color','k','linewidth',0.1);
                hold off
                % Annotate plot
                title(sprintf('%s Bin%i',ID.Circuit{Idx_Circuit},Idx_Bin));
                xlabel(ID.Label{2});
                ylabel(ID.Label{Idx_Color});
                ax = gca;
                ax.XScale = 'log'; xlim([1 2e5]);
                ax.YScale = 'log'; ylim([1 2e5]);
                grid on
            end
            
            % Make plot screen size and save it as .png
            fig.Units = 'normalized';
            fig.OuterPosition = [0 0 1 1];
            fig.PaperPositionMode = 'auto';
            fig.Visible = 'on';
%           % Save figure as .png-file
            print(fullfile(folder_bubble,sprintf('Input_%s_Bin%i_Mean.png',ID.Label{Idx_Color},Idx_Bin)),'-dpng');
            % Save figure as .fig-file
            savefig(fullfile(folder_bubble,sprintf('Input_%s_Bin%i_Mean.fig',ID.Label{Idx_Color},Idx_Bin)));
            close all hidden;
        end
    end
end

%% Composite plots

for Idx_Rep = 1:ID.Replicates
    % Output folder
    folder_composite = fullfile(folder,ID.ReplicateField{Idx_Rep},'Composite_Plots');
    mkdir(folder_composite);
    
    for Idx_Bin = 1:1:ID.TMBins
        fig = figure;
        for Idx_Out = 3:ID.NrOut+2
            for Idx_Circuit = 1:ID.Circuits
                subplot(double(ceil(sqrt(ID.NrOut*ID.Circuits))),double(ceil(sqrt(ID.NrOut*ID.Circuits))),double(Idx_Circuit+(Idx_Out-3)*ID.Circuits))
                % Create flow plots
                cloudPlot(Data.(ID.ReplicateField{Idx_Rep}).comp_biex{Idx_Circuit,Idx_Bin}(:,2),Data.(ID.ReplicateField{Idx_Rep}).comp_biex{Idx_Circuit,Idx_Bin}(:,Idx_Out),[0 4.5 0 4.5],false,[250,250]);
                hold on
                % Add the PFAFF output ontop
                Subidx_Rep = [[1:ID.InputLvl]+(Idx_Rep-1)*ID.InputLvl,...
                    [ID.InputLvl*ID.Replicates*1+1:ID.InputLvl*ID.Replicates*1+ID.InputLvl]+(Idx_Rep-1)*ID.InputLvl,...
                    [ID.InputLvl*ID.Replicates*2+1:ID.InputLvl*ID.Replicates*2+ID.InputLvl]+(Idx_Rep-1)*ID.InputLvl,...
                    [ID.InputLvl*ID.Replicates*3+1:ID.InputLvl*ID.Replicates*3+ID.InputLvl]+(Idx_Rep-1)*ID.InputLvl];
                scatter(BiEx(DataSet.X{Idx_Circuit,Idx_Bin}(Subidx_Rep),LookUpTable,PlotRange),BiEx(DataSet.Y{Idx_Out}{Idx_Circuit,Idx_Bin}(Subidx_Rep),LookUpTable,PlotRange),100*DataSet.W{Idx_Out}{Idx_Circuit,Idx_Bin}(Subidx_Rep),'markeredgecolor','k');
                hold off
                % Use Material Design colormap
                colormap(MaterialMap); caxis([0,100]);
                % Annotate plot
                title(sprintf('%s, Bin %i',ID.Circuit{Idx_Circuit},Idx_Bin));
                xlabel(ID.Label{2});
                ylabel(ID.Label{Idx_Out});
                ax = gca;
                ax.XTick = ID.Flow_scale; ax.XTickLabel = {'0','10^2','10^3','10^4','10^5'};
                ax.YTick = ID.Flow_scale; ax.YTickLabel = {'0','10^2','10^3','10^4','10^5'};
                grid on
            end
        end
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1.08/1.92 1];
        fig.PaperPositionMode = 'auto';
        fig.Visible = 'on';
%       % Save figure as .png-file
        print(fullfile(folder_composite,sprintf('CompositeFlow_Plots_Rep%i_Bin%i.png',Idx_Rep,Idx_Bin)),'-dpng');
        % Save figure as .fig-file
        savefig(fullfile(folder_composite,sprintf('CompositeFlow_Plots_Rep%i_Bin%i.fig',Idx_Rep,Idx_Bin)));
        close all hidden;
    end
end

%% Save workspace and variables
disp('Analysis finished. Saving ...');
save(fullfile(folder,'DataSet.mat'),'DataSet');
save(fullfile(folder,'Analysis.mat'),'Data','DataSet','ID','-v7.3');

%% Stop timer *********************************************************** %
disp(['Total duration: ',num2str(toc/60),' minutes']);
end

%% Additional functions
function [binned_events,PeakValues] = PeakBinning(array,flagBin,Idx_Input,InPeakValue,LookUpTable,PlotRange)
global BinEdge BinCenter

switch flagBin
    case 1
        % Case 1: classic event single binning
        Bin_number = 10;
        
        % Sort descending in the transfection control column
        s_array = sortrows(array,-1);
        
        % Set intensity range
        lb = prctile(s_array(:,1),2.5);
        ub = prctile(s_array(:,1),97.5);
        
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
        
        % Dummy variable
        PeakValues = NaN;
        
    case 2
        % Case 2: Bin according to the peak of the TM distribution
        % Transform data into biex space
        biex_array = BiEx(array(:,1),LookUpTable,PlotRange);
        
        % Generate a histogram of the TM data
        N = histcounts(biex_array,BinEdge);
        
        % Fitoptions
        opt1 = fitoptions('gauss1','Lower',[0 0 0],'Upper',[Inf 4.5 4.5]); %,'StartPoint',[max(N),2.1,1]);
        opt2 = fitoptions('gauss2','Lower',[0 0 0 0 0 0],'Upper',[Inf 4.5 4.5 Inf 4.5 4.5]); %,'StartPoint',[max(N),0.6,1,max(N),2.1,1]);
        
        % Fit Gauss1 to data
        [f1,gof1] = fit(BinCenter,N','gauss1', opt1);
        % Fit gauss2 to data
        [f2,~] = fit(BinCenter,N','gauss2', opt2);
        if gof1.rsquare > 0.95
            Mode = f1.b1;
        else
            Mode = max([f2.b1,f2.b2]);
        end
        
        % Sort biex array
        sort_array = sortrows(biex_array);
        % Determine the percentile of the peak
        [~,i] = min(abs(sort_array-Mode));
        PeakValues = i/length(sort_array);
        
        % Set boundaries
        HalfWindow = 0.075; % Half window size
        lb = PeakValues - HalfWindow;
        ub = PeakValues + HalfWindow;
        
        % Select events
        Idx_Bin = biex_array >= prctile(biex_array,lb*100) & ...
            biex_array <= prctile(biex_array,ub*100);
        binned_events{1} = array(Idx_Bin,:);
        
    case 3
        % Case 3: Bin according to the percentile range of the Peak from 0 Input
        if Idx_Input == 1
            % Transform data into biex space
            biex_array = BiEx(array(:,1),LookUpTable,PlotRange);
            
            % Generate a histogram of the TM data
            N = histcounts(biex_array,BinEdge);
            
            % Fitoptions
            opt1 = fitoptions('gauss1','Lower',[0 0 0],'Upper',[Inf 4.5 4.5]); %,'StartPoint',[max(N),2.1,1]);
            opt2 = fitoptions('gauss2','Lower',[0 0 0 0 0 0],'Upper',[Inf 4.5 4.5 Inf 4.5 4.5]); %,'StartPoint',[max(N),0.6,1,max(N),2.1,1]);
            
            % Fit Gauss1 to data
            [f1,gof1] = fit(BinCenter,N','gauss1', opt1);
            % Fit gauss2 to data
            [f2,~] = fit(BinCenter,N','gauss2', opt2);
            if gof1.rsquare > 0.95
                Mode = f1.b1;
            else
                Mode = max([f2.b1,f2.b2]);
            end
            
            % Sort biex array
            sort_array = sortrows(biex_array);
            % Determine the percentile of the peak
            [~,i] = min(abs(sort_array-Mode));
            PeakValues = i/length(sort_array);
            
            % Set boundaries
            HalfWindow = 0.075; % Half window size
            lb = PeakValues - HalfWindow;
            ub = PeakValues + HalfWindow;
            
            % Select events
            Idx_Bin = biex_array >= prctile(biex_array,lb*100) & ...
                biex_array <= prctile(biex_array,ub*100);
            binned_events{1} = array(Idx_Bin,:);
        else
            % Use PeakValue from input sample #1
            PeakValues = InPeakValue(1);
            
            % Transform data into biex space
            biex_array = BiEx(array(:,1),LookUpTable,PlotRange);
            
            % Set boundaries
            HalfWindow = 0.075; % Half window size
            lb = PeakValues - HalfWindow;
            ub = PeakValues + HalfWindow;
            
            % Select events
            Idx_Bin = biex_array >= prctile(biex_array,lb*100) & ...
                biex_array <= prctile(biex_array,ub*100);
            binned_events{1} = array(Idx_Bin,:);
        end
        
    case 4
        % Case 4: classic event single binning for stable
        Bin_number = 1;
        
        % Sort descending in the transfection control column
        s_array = sortrows(array,-1);
        
        % Set intensity range
        lb = prctile(s_array(:,1),2.5);
        ub = prctile(s_array(:,1),97.5);
        
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
        
        % Dummy variable
        PeakValues = NaN;
end
end