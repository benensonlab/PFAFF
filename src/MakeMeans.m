function [DataSet] = MakeMeans(DataSet,ID)
% Uses the DataSet dollection of variables and averages them in biex space.
% The output of the function is DataSet_bar that includes the nanmean of X,
% Y and A, as well as the SD (lower and upper bounds) and a new weight
% vector W.

load LookUpTable

% Determine peak vector (peaks x replicates)
Size_PR = [1:ID.Replicates;ID.Replicates+1:ID.Replicates*2;2*ID.Replicates+1:ID.Replicates*3;3*ID.Replicates+1:ID.Replicates*4];

for Idx_Circuit = 1:ID.Circuits
    for Idx_Bin = 1:ID.TMBins
        % Reshape X vector
        X{Idx_Circuit,Idx_Bin} = reshape(BiEx(DataSet.X{Idx_Circuit,Idx_Bin},LookUpTable,PlotRange),ID.InputLvl,4*ID.Replicates); % 12 = 4 peaks*3 Reps
        % Compute nanmean for X values
        X_bar{Idx_Circuit,Idx_Bin}(:,1) = nanmean(X{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),2);
        X_bar{Idx_Circuit,Idx_Bin}(:,2) = nanmean(X{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),2);
        X_bar{Idx_Circuit,Idx_Bin}(:,3) = nanmean(X{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),2);
        X_bar{Idx_Circuit,Idx_Bin}(:,4) = nanmean(X{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),2);
        % Compute nanstd for X values
        X_SD{Idx_Circuit,Idx_Bin}(:,1) = nanstd(X{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),[],2);
        X_SD{Idx_Circuit,Idx_Bin}(:,2) = nanstd(X{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),[],2);
        X_SD{Idx_Circuit,Idx_Bin}(:,3) = nanstd(X{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),[],2);
        X_SD{Idx_Circuit,Idx_Bin}(:,4) = nanstd(X{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),[],2);
        % Determine lower and upper bounds of SD and retransform
        X_SD_lb{Idx_Circuit,Idx_Bin} = reshape(invBiEx(X_bar{Idx_Circuit,Idx_Bin}-X_SD{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
        X_SD_ub{Idx_Circuit,Idx_Bin} = reshape(invBiEx(X_bar{Idx_Circuit,Idx_Bin}+X_SD{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
        % Reshape into original layout and retransform
        X_bar{Idx_Circuit,Idx_Bin} = reshape(invBiEx(X_bar{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
        
        for Idx_Color = 3:ID.NrOut+2
            % Reshape Y vector
            Y{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(BiEx(DataSet.Y{Idx_Color}{Idx_Circuit,Idx_Bin},LookUpTable,PlotRange),ID.InputLvl,4*ID.Replicates);
            % Compute nanmean for Y values
            Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1) = nanmean(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),2);
            Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2) = nanmean(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),2);
            Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,3) = nanmean(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),2);
            Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,4) = nanmean(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),2);
            % Compute nanstd for Y values
            Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1) = nanstd(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),[],2);
            Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2) = nanstd(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),[],2);
            Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}(:,3) = nanstd(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),[],2);
            Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}(:,4) = nanstd(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),[],2);
            % Determine lower and upper bounds of SD and retransform
            Y_SD_lb{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(invBiEx(Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}-Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
            Y_SD_ub{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(invBiEx(Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}+Y_SD{Idx_Color}{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
            % Reshape into original layout and retransform
            Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(invBiEx(Y_bar{Idx_Color}{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
            
            % Generate new W_bar vector
            W{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(BiEx(DataSet.W{Idx_Color}{Idx_Circuit,Idx_Bin},LookUpTable,PlotRange),ID.InputLvl,4*ID.Replicates);
            % Compute nanmean for Y values
            W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1) = nanmean(W{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),2);
            W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2) = nanmean(W{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),2);
            W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,3) = nanmean(W{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),2);
            W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,4) = nanmean(W{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),2);
            % Reshape into original layout and retransform
            W_bar{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(invBiEx(W_bar{Idx_Color}{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
            
            % % Reshape A vector
            % A{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(BiEx(DataSet.A{Idx_Color}{Idx_Circuit,Idx_Bin},LookUpTable,PlotRange),ID.InputLvl,4*ID.Replicates);
            % % Compute nanmean for A values
            % A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1) = nanmean(A{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(1,:)),2);
            % A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2) = nanmean(A{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(2,:)),2);
            % A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,3) = nanmean(A{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(3,:)),2);
            % A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,4) = nanmean(A{Idx_Color}{Idx_Circuit,Idx_Bin}(:,Size_PR(4,:)),2);
            % % Reshape into original layout and retransform
            % A_bar{Idx_Color}{Idx_Circuit,Idx_Bin} = reshape(invBiEx(A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}),2*ID.InputLvl,2);
            %
            % % Generate new W_bar vector
            % W_bar{Idx_Color}{Idx_Circuit,Idx_Bin} = [A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1)./nansum(A_bar{Idx_Color}{Idx_Circuit,Idx_Bin},2),...
            %   A_bar{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2)./nansum(A_bar{Idx_Color}{Idx_Circuit,Idx_Bin},2)];
        end
    end
end

% Write into structure
% DataSet.A_bar = A_bar;

DataSet.X_bar = X_bar;
DataSet.X_SD_lb = X_SD_lb;
DataSet.X_SD_ub = X_SD_ub;

DataSet.Y_bar = Y_bar;
DataSet.Y_SD_lb = Y_SD_lb;
DataSet.Y_SD_ub = Y_SD_ub;

DataSet.W_bar = W_bar;