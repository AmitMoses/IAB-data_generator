classdef IABdonor < CommUnit

    properties
        UE_con  % list of connected UE
        UE_con_rssi
        UE_con_CQI
        BW
        MaxUEperBS = 20;
        Tx = WaveForm
        ResourceGrid_DL
        ResourceGrid_UL
    end
    
    methods
        %%
        function obj = Connect2UE(obj,UE_list, Pathloss,CQI2SNR)
%             disp('recived')
            % Calculate RSSI for every link to UE
            Psens = -Pathloss(obj.UE_con) + 3 + 3 +  [UE_list(obj.UE_con).Ptx]; % Calculate reviced power from each user
            obj.UE_con_rssi = Psens;
            
            % Calculate CQI for every link to UE
            if isempty(Psens) == 0              
                N0 = -174 + 10*log10(obj.BW);
                SNR = Psens - N0;
                CQI_idx = size(CQI2SNR,1) + 1 - sum(double(CQI2SNR(:,2)>SNR),1) - 2;
                obj.UE_con_CQI = CQI2SNR(CQI_idx+1,1).';
%                 obj.UE_con_CQI = CQI_idx;
            end
        end
        
        %%
        function obj = resource_allocation(obj)
            obj.Tx = setup(obj.Tx, 10e6, 2);
%             RB_Matrix = zeros(obj.Tx.NumPRB, obj.Tx.NumSlotPerSubframe);
            strategy = 'Random';
            switch strategy
                case 'Random'
                   DL_RB_Matrix_index = randi(length(obj.UE_con), obj.Tx.NumPRB,obj.Tx.NumSlotPerSubframe ) ;
                   obj.ResourceGrid_DL = obj.UE_con(DL_RB_Matrix_index);
                   UL_RB_Matrix_index = randi(length(obj.UE_con), obj.Tx.NumPRB,obj.Tx.NumSlotPerSubframe ) ;
                   obj.ResourceGrid_UL = obj.UE_con(UL_RB_Matrix_index);
            end
        end
        
    end
    
end