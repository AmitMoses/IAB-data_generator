classdef UE < CommUnit
    % User Equpment
    properties
        BS_con_id       % Servcie providing BS ID
        BS_con_rssi     % Servcie providing BS RSSI
        BS_con_CQI      % Servcie providing BS CQI
        BW
        Tx = zeros(24,14);
        hops
    end
    
    methods 
        %%
        function [net] = Connect2BS(obj, net, BS_pathloss, IAB_pathloss, IABidcator)
            % UE search best BS to Connect
            global IABnode_num
            global CQI2SNR
            
            if 1    % make each IAB to be with max of 'MaxUEperBS' UE 
                if length(net.IABdonors(:).UE_con) == net.IABdonors.MaxUEperBS && not(IABidcator)
                    BS_pathloss = inf;
                end
                for i = 1:IABnode_num
                   if length(net.IABnodes(i).gNB.UE_con) == net.IABnodes(i).gNB.MaxUEperBS
                       IAB_pathloss(i) = inf;
                   end
                end
            end
            
            % choose between IAB-Donor to IAB-node
            if IABidcator == 1
                [rssi, win_idx, Psens_BS, Psens_IAB, max_cqi_idx] = chooss_gNB_cqi(obj, net, BS_pathloss, IAB_pathloss, IABidcator);
            else
                [rssi, win_idx, Psens_BS, Psens_IAB, max_cqi_idx] = choose_gNb(obj, net, BS_pathloss, IAB_pathloss, IABidcator);
            end
            
            
                       
            % CQI Calculation
            CQI = rssi2cqi(obj,rssi);
            
%             CQI = find(CQIvsRSSI(:,2)<rssi,1);
            if CQI == CQI2SNR(1,1)
                % UE unsuccsses to connect
                disp('not connected')
                net.disconnected_users = [net.disconnected_users obj.ID];
            else
                % UE succsses to connect
                if win_idx == 1
%                     [obj.BS_con_rssi,BS_con_idx] = max(Psens_BS);           % find best Node to connect;
                    obj.BS_con_rssi = rssi;
                    BS_con_idx = max_cqi_idx;
                    obj.BS_con_id = net.IABdonors(BS_con_idx).ID;
                    obj.BS_con_CQI = CQI;
                    net.IABdonors(BS_con_idx).UE_con = [ net.IABdonors(BS_con_idx).UE_con obj.ID ];  % Tell the BS that connection is established
                % UE succsses to connect
                elseif win_idx == 2
%                     [obj.BS_con_rssi,BS_con_idx] = max(Psens_IAB);           % find best Node to connect;
                    obj.BS_con_rssi = rssi;
                    BS_con_idx = max_cqi_idx;
                    obj.BS_con_id = net.IABnodes(BS_con_idx).gNB.ID;
                    obj.BS_con_CQI = CQI;
                    net.IABnodes(BS_con_idx).gNB.UE_con = [ net.IABnodes(BS_con_idx).gNB.UE_con obj.ID ];  % Tell the BS that connection is established
                end
            end
            % insert change of UE to net object
            if IABidcator == 0
                net.users(obj.ID) = obj;
            else
                net.IABnodes(obj.ID - length(net.users)).UE = obj;
                net.IABnodes(obj.ID - length(net.users)).connect2Donor = 1;
            end
        end
        
        %% choose_gNb base on rssi
        function [rssi, win_idx, Psens_BS, Psens_IAB, max_cqi_idx] = choose_gNb(obj, net, BS_pathloss, IAB_pathloss, IABidcator)
            global CQI2SNR
            % Recive power from IAB-Donor
            Psens_BS = -BS_pathloss + 3 + 3 +  [net.IABdonors(:).Ptx].';
            [rssi_BS,rssi_BS_idx] = max(Psens_BS);           % find best Node to connect;
            
            % Recive power from IAB-Nodes
            gNB_list = [net.IABnodes(:).gNB];
            Psens_IAB = -IAB_pathloss + 3 + 3 +  [gNB_list(:).Ptx].';
            if IABidcator
                is_cottected = [net.IABnodes(:).connect2Donor];
                Psens_IAB(is_cottected==0)=-inf;
            end
            [rssi_IAB,rssi_IAB_idx] = max(Psens_IAB);
            
            % choose between IAB-Donor to IAB-node
            [rssi,win_idx] = max([rssi_BS,rssi_IAB]);
            if win_idx == 1
                max_cqi_idx = rssi_BS_idx;
            elseif win_idx == 2
                max_cqi_idx = rssi_IAB_idx;
            end
            
        end
        
        %% choose_gNb base on CQI
        function [rssi, win_idx, Psens_BS, Psens_IAB, max_cqi_idx] = chooss_gNB_cqi(obj, net, BS_pathloss, IAB_pathloss, IABidcator)
            % Recive power from IAB-Donor
            Psens_BS = -BS_pathloss + 3 + 3 +  [net.IABdonors(:).Ptx].';
            [rssi_BS,~] = max(Psens_BS);           % find best Node to connect;

            % Recive power from IAB-Nodes
            gNB_list = [net.IABnodes(:).gNB];
            Psens_IAB = -IAB_pathloss + 3 + 3 +  [gNB_list(:).Ptx].';
            if IABidcator
                is_cottected = [net.IABnodes(:).connect2Donor];
                Psens_IAB(is_cottected==0)=-inf;
            end
            [rssi_IAB,~] = max(Psens_IAB);

            % choose between IAB-Donor to IAB-node
            CQI_BS = rssi2cqi(obj, Psens_BS);
            CQI_IAB = rssi2cqi(obj, Psens_IAB.');
            concat_CQI = [CQI_BS,CQI_IAB];
            [~,max_cqi_idx] = max(concat_CQI);
            if max_cqi_idx == 1
                rssi = Psens_BS(max_cqi_idx);
                win_idx = 1;
            else
                max_cqi_idx = max_cqi_idx-length(Psens_BS);
                rssi = Psens_IAB(max_cqi_idx);
                win_idx  = 2;
            end
              
        end
        
        %% RSSI to CQI
        function [CQI] = rssi2cqi(obj, rssi)
            global CQI2SNR
            N0 = -174 + 10*log10(obj.BW);
            SNR = rssi - N0;
%             CQI = find(CQI2SNR(:,2)>SNR,1) - 2;
            
            map = bsxfun(@gt,CQI2SNR(:,2),SNR);
            [~, firstoverthresh]=max(map,[],1);
            CQI = firstoverthresh - 2;
        end
                       
    end
end
              
