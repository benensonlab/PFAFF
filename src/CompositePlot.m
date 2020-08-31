function CompositePlot(Data,DataSet,ID)

% Labels
DOX = [0,1,3.5,12,42,144,500,1500];
circ = {'A','B','C','D','E','F','G','H'};
FP = {'BFP','Cerulean','Citrine','Cherry'};
cm = generateColormap;
load LookUpTable
FACS_scale = BiEx([0,100,1000,10000,100000],LookUpTable,PlotRange);
Repfield = fieldnames(Data);

for REP = 1:3
    for Circuit = 1:8
        for Bin = 1:10
            for Dox = 1:8
                ds{Circuit,Bin} = vertcat(Data.(Repfield{REP}).biex_data{Circuit,:,Bin});
            end;
        end;
    end;
        
    folder = fullfile(date,'Data');
    mkdir(folder);
    
    for Bin = 1:10
        figure
        for Circuit = 1:8
            subplot(4,4,Circuit)
            cloudPlot(ds{Circuit,Bin}(:,2),ds{Circuit,Bin}(:,3),[0 4.5 0 4.5],false,[250,250]);
            hold on
            Idx_Rep = [[1:8]+REP*8-8,[25:32]+REP*8-8,[49:56]+REP*8-8,[73:80]+REP*8-8];
            scatter(BiEx(DataSet.X{Circuit,Bin}(Idx_Rep),LookUpTable,PlotRange),BiEx(DataSet.Y{3}{Circuit,Bin}(Idx_Rep),LookUpTable,PlotRange),100*DataSet.W{3}{Circuit,Bin}(Idx_Rep),'markeredgecolor','k');
            hold off
            colormap(cm); caxis([0,100]);
            xlabel(FP{4});
            ylabel(FP{2});
            title(sprintf('Circuit %c, Bin %i',circ{Circuit},Bin));
            ax = gca;
            ax.XTick = FACS_scale; ax.XTickLabel = {'0','10^2','10^3','10^4','10^5'};
            ax.YTick = FACS_scale; ax.YTickLabel = {'0','10^2','10^3','10^4','10^5'};
            grid on
            
            subplot(4,4,Circuit+8)
            cloudPlot(ds{Circuit,Bin}(:,2),ds{Circuit,Bin}(:,4),[0 4.5 0 4.5],false,[250,250]);
            hold on
            scatter(BiEx(DataSet.X{Circuit,Bin}(Idx_Rep),LookUpTable,PlotRange),BiEx(DataSet.Y{4}{Circuit,Bin}(Idx_Rep),LookUpTable,PlotRange),100*DataSet.W{4}{Circuit,Bin}(Idx_Rep),'markeredgecolor','k');
            hold off
            colormap(cm); caxis([0,100]);
            xlabel(FP{4});
            ylabel(FP{3});
            title(sprintf('Circuit %c, Bin %i',circ{Circuit},Bin));
            ax = gca;
            ax.XTick = FACS_scale; ax.XTickLabel = {'0','10^2','10^3','10^4','10^5'};
            ax.YTick = FACS_scale; ax.YTickLabel = {'0','10^2','10^3','10^4','10^5'};
            grid on
        end;
        fig = gcf;
        fig.Units = 'normalized';fig.OuterPosition = [0 0 1.08/1.92 1];fig.PaperPositionMode = 'auto';
        print(fullfile(folder,sprintf('CompositeFACS_Plots_Data_Rep%i_Bin%i.png',REP,Bin)),'-dpng','-r0');
    end;
    close all hidden
end;

