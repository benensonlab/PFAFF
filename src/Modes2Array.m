function [DataSet] = Modes2Array(data,ID)
% Extracts modes from the gaussian fit into single variables (a: height, b:
% mode position, c: width). Variables are 5D cell arrays, and indices are:
% e.g. b{Color,InputPeak,OutputPeak,Circuit,Bin};

for Idx_Rep = 1:ID.Replicates
    for Idx_Circuit = 1:ID.Circuits
        for Idx_Input = 1:ID.InputLvl
            for Idx_Bin = 1:ID.TMBins
                for InPeak = 1:2
                    % Input 
                    if length(data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(:,InPeak)) < 3
                        a{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(1,InPeak);
                        b{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,InPeak));
                        c{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = NaN;
                        area{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = NaN;
                    else
                        a{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(1,InPeak);
                        b{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(2,InPeak));
                        c{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(3,InPeak));
                        area{2,InPeak,1,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Input{Idx_Circuit,Idx_Input,Idx_Bin}(4,InPeak);
                    end
                    
                    % Output
                    for Idx_Out = 1:ID.NrOut
                        for OutPeak = 1:2
                            if length(data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(:,InPeak)) < 3
                                a{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(1,OutPeak);
                                b{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(2,OutPeak));
                                c{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = NaN;
                                area{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = NaN;
                            else
                                a{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(1,OutPeak);
                                b{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(2,OutPeak));
                                c{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = invBiEx(data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(3,OutPeak));
                                area{Idx_Out+2,InPeak,OutPeak,Idx_Circuit,Idx_Bin}(Idx_Input+(Idx_Rep-1)*ID.InputLvl,:) = data.(ID.ReplicateField{Idx_Rep}).Mode_Output{Idx_Out}{InPeak,Idx_Circuit,Idx_Input,Idx_Bin}(4,OutPeak);
                            end
                        end
                    end
                end
            end
        end
    end
end


for Idx_Circuit = 1:ID.Circuits
    for Idx_Bin = 1:ID.TMBins
        for Idx_Color = 3:ID.NrOut+2
            A{Idx_Color}{Idx_Circuit,Idx_Bin} = [a{Idx_Color,1,1,Idx_Circuit,Idx_Bin},a{Idx_Color,1,2,Idx_Circuit,Idx_Bin};a{Idx_Color,2,1,Idx_Circuit,Idx_Bin},a{Idx_Color,2,2,Idx_Circuit,Idx_Bin}];
            Area{Idx_Color}{Idx_Circuit,Idx_Bin} = [area{Idx_Color,1,1,Idx_Circuit,Idx_Bin},area{Idx_Color,1,2,Idx_Circuit,Idx_Bin};area{Idx_Color,2,1,Idx_Circuit,Idx_Bin},area{Idx_Color,2,2,Idx_Circuit,Idx_Bin}];
            Y{Idx_Color}{Idx_Circuit,Idx_Bin} = [b{Idx_Color,1,1,Idx_Circuit,Idx_Bin},b{Idx_Color,1,2,Idx_Circuit,Idx_Bin};b{Idx_Color,2,1,Idx_Circuit,Idx_Bin},b{Idx_Color,2,2,Idx_Circuit,Idx_Bin}];
            W{Idx_Color}{Idx_Circuit,Idx_Bin} = [Area{Idx_Color}{Idx_Circuit,Idx_Bin}(:,1)./nansum(Area{Idx_Color}{Idx_Circuit,Idx_Bin},2),Area{Idx_Color}{Idx_Circuit,Idx_Bin}(:,2)./nansum(Area{Idx_Color}{Idx_Circuit,Idx_Bin},2)];
            
            for i = 1:length(Y{Idx_Color}{Idx_Circuit,Idx_Bin}(:))/2
                if W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) >= W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2)
                    S{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) = Y{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1);
                    
                elseif W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) < W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2)
                    S{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) = Y{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2);
                    
                elseif isnan(W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1)) && ~isnan(W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2))
                    S{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) = Y{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2);
                    
                elseif ~isnan(W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1)) && isnan(W{Idx_Color}{Idx_Circuit,Idx_Bin}(i,2))
                    S{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) = Y{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1);
                    
                else
                    S{Idx_Color}{Idx_Circuit,Idx_Bin}(i,1) = NaN;
                end
            end
            
        end
        X{Idx_Circuit,Idx_Bin} = [b{2,1,1,Idx_Circuit,Idx_Bin},b{2,1,1,Idx_Circuit,Idx_Bin};b{2,2,1,Idx_Circuit,Idx_Bin},b{2,2,1,Idx_Circuit,Idx_Bin}];
    end
end

% DataSet.a = a;
% DataSet.b = b;
% DataSet.c = c;
% DataSet.area = area;

% DataSet.A = A;
% DataSet.Area = Area;
% DataSet.S = S;
DataSet.W = W;
DataSet.X = X;
DataSet.Y = Y;