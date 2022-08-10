classdef saving
    properties
        ue_data
        iab_data
    end
    methods
        %% Setup object properties
        function obj = set(obj, UE_database, IAB_database)
            obj.ue_data = UE_database;
            obj.iab_data = IAB_database;
        end
        
        %% update one row of ue_data 
        function obj = update_ue(obj, net, row)
            global Ue_Num
            for i=0:Ue_Num-1
                idxDL = find(net.Topology.Edges.EndNodes(:,1) == i+1);
                idxUL = find(net.Topology.Edges.EndNodes(:,2) == i+1);
                obj.ue_data(row,i*5+1) = net.users(i+1).BS_con_id;               % connected BS
                obj.ue_data(row,i*5+2) = net.Topology.Edges.Capacity(idxDL);     % DL_app
                obj.ue_data(row,i*5+3) = net.Topology.Edges.CQI(idxDL);          % DL_CQI
                obj.ue_data(row,i*5+4) = net.Topology.Edges.Capacity(idxUL);     % UL_app
                obj.ue_data(row,i*5+5) = net.Topology.Edges.CQI(idxUL);          % UL_CQI
            end 
        end
        
        %% Updata one row of iab_data
        function obj = update_iab(obj, net, row)
            
            global Ue_Num
            global IABnode_num
            global IABdonor_Num
            global max_bachaul_num
            
            for i=0:IABnode_num -1
                 IABnode_ID = net.IABnodes(i+1).ID;
                 connected_gNB = net.IABnodes(i+1).UE.BS_con_id;

                 for backhaul_idx = 0:length(connected_gNB)-1
                     index = 5*(i*max_bachaul_num + backhaul_idx);
                     backhaul = connected_gNB(backhaul_idx + 1);
                     idxDL = find(ismember(net.Topology.Edges.EndNodes, [backhaul, IABnode_ID],'rows'));
                     idxUL = find(ismember(net.Topology.Edges.EndNodes, [IABnode_ID, backhaul],'rows'));

                    if isnan(net.IABnodes(i+1).UE.BS_con_id) 
                        obj.iab_data(row,index+1) = 0;
                    else
                        obj.iab_data(row,index+1) = net.IABnodes(i+1).UE.BS_con_id(backhaul_idx + 1);          % Connected BS
                    end

                    if isnan(net.Topology.Edges.Capacity(idxDL))
                        obj.iab_data(row,index+2) = 0;
                    else
                        obj.iab_data(row,index+2) = net.Topology.Edges.Capacity(idxDL);      % DL_app
                    end

                    if isnan(net.Topology.Edges.CQI(idxDL))
                        obj.iab_data(row,index+3) = 0;
                    else
                        obj.iab_data(row,index+3) = net.Topology.Edges.CQI(idxDL);           % DL_CQI
                    end

                    if isnan(net.Topology.Edges.Capacity(idxUL))
                        obj.iab_data(row,index+4) = 0;
                    else
                        obj.iab_data(row,index+4) = net.Topology.Edges.Capacity(idxUL);      % UL_app
                    end

                    if isnan(net.Topology.Edges.CQI(idxUL))
                        obj.iab_data(row,index+5) = 0;
                    else
                        obj.iab_data(row,index+5) = net.Topology.Edges.CQI(idxUL);           % UL_CQI
                    end

                 end
            end

            % saving IAB-Donors database into IAB_database matrix
            for i=0:IABdonor_Num -1
                index = (i + IABnode_num*max_bachaul_num)*5;
                 idxDL = intersect(...
                    find(net.Topology.Edges.EndNodes(:,2) <= Ue_Num ),...
                    find(net.Topology.Edges.EndNodes(:,1) == net.IABdonors(i+1).ID));
                 idxUL = intersect(...
                    find(net.Topology.Edges.EndNodes(:,1) <= Ue_Num ),...
                    find(net.Topology.Edges.EndNodes(:,2) == net.IABdonors(i+1).ID));
                obj.iab_data(row,index+1) = -1;          % Connected BS
                if isnan(sum(net.Topology.Edges.Capacity(idxDL))) 
                    obj.iab_data(row,index+2) = 0;
                else
                    obj.iab_data(row,index+2) = sum(net.Topology.Edges.Capacity(idxDL));         % DL_app
                end

                if isnan(floor(mean(net.Topology.Edges.CQI(idxDL))))
                    obj.iab_data(row,index+3) = 0;
                else
                    obj.iab_data(row,index+3) = floor(mean(net.Topology.Edges.CQI(idxDL)));      % DL_CQI
                end

                if isnan(sum(net.Topology.Edges.Capacity(idxUL)))
                    obj.iab_data(row,index+4) = 0;
                else
                    obj.iab_data(row,index+4) = sum(net.Topology.Edges.Capacity(idxUL));         % UL_app
                end

                if isnan(floor(mean(net.Topology.Edges.CQI(idxUL))))
                    obj.iab_data(row,index+5) = 0;
                else
                    obj.iab_data(row,index+5) = floor(mean(net.Topology.Edges.CQI(idxUL)));      % UL_CQI
                end
            end 
        end
        
        %% Save the UE data
        function UE_Table(obj)
            global Ue_Num
            
            % Generate UE database labels
            UE_database_labels = [];
            for  i=1:Ue_Num
                str_0 = ['UE', num2str(i),'-Con_BS '];
                str_1 = ['UE', num2str(i),'-DL_app '];
                str_2 = ['UE', num2str(i),'-DL_CQI '];
                str_3 = ['UE', num2str(i),'-UL_app '];
                str_4 = ['UE', num2str(i),'-UL_CQI '];
                UE_database_labels = [UE_database_labels,str_0, str_1,str_2,str_3,str_4];
            end
            UE_database_labels = split(UE_database_labels);
            UE_database_labels = UE_database_labels(1:end-1).';

            % Convert database matrix to table
            UE_database_Table = table(obj.ue_data);
            UE_database_Table = splitvars(UE_database_Table);
            UE_database_Table.Properties.VariableNames = UE_database_labels;
            
            disp('Saving UE data into csv file...')
            writetable(UE_database_Table,'UE_database.csv')  
            disp('Done.')
        end
        
        %% Save the IAB data
        function IAB_Table(obj)
            global Ue_Num
            global IABnode_num
            global IABdonor_Num
            global max_bachaul_num
            
            % Generate IAB database labels
            IAB_database_labels = [];
                for  i=1:IABnode_num
                    for j=1:max_bachaul_num
                        str_0 = ['IAB', num2str(i+Ue_Num),['-Con_BS_',num2str(j)],' '];
                        str_1 = ['IAB', num2str(i+Ue_Num),['-DL_app_',num2str(j)],' '];
                        str_2 = ['IAB', num2str(i+Ue_Num),['-DL_CQI_',num2str(j)],' '];
                        str_3 = ['IAB', num2str(i+Ue_Num),['-UL_app_',num2str(j)],' '];
                        str_4 = ['IAB', num2str(i+Ue_Num),['-UL_CQI_',num2str(j)],' '];
                        IAB_database_labels = [IAB_database_labels,str_0, str_1,str_2,str_3,str_4];
                    end
                end

                for  i=1+IABnode_num:IABnode_num + IABdonor_Num
                    str_0 = ['IAB', num2str(i+Ue_Num),'-Con_BS '];
                    str_1 = ['IAB', num2str(i+Ue_Num),'-DL_app '];
                    str_2 = ['IAB', num2str(i+Ue_Num),'-DL_CQI '];
                    str_3 = ['IAB', num2str(i+Ue_Num),'-UL_app '];
                    str_4 = ['IAB', num2str(i+Ue_Num),'-UL_CQI '];
                    IAB_database_labels = [IAB_database_labels,str_0, str_1,str_2,str_3,str_4];
                end
                IAB_database_labels = split(IAB_database_labels);
                IAB_database_labels = IAB_database_labels(1:end-1).';

                % Convert database matrix to table
                IAB_database_Table = table(obj.iab_data);
                IAB_database_Table = splitvars(IAB_database_Table);
                IAB_database_Table.Properties.VariableNames = IAB_database_labels;
                
                disp('Saving IAB data into csv file...')
                writetable(IAB_database_Table,'IAB_database.csv')
                disp('Done.')
        end
        
        %% Save the IAB_graph data
        function IABgraph(obj, net_obj, iter)
           
            global Ue_Num
            global IABnode_num
            global IABdonor_Num
            
            idxDL = intersect(...
                find(net_obj.Topology.Edges.EndNodes(:,2) <= Ue_Num ),...
                find(net_obj.Topology.Edges.EndNodes(:,1) == net_obj.IABdonors(IABdonor_Num).ID));
            idxUL = intersect(...
                find(net_obj.Topology.Edges.EndNodes(:,1) <= Ue_Num ),...
                find(net_obj.Topology.Edges.EndNodes(:,2) == net_obj.IABdonors(IABdonor_Num).ID));
            
            idx =(1:IABnode_num+IABdonor_Num)+Ue_Num;
            H = subgraph(net_obj.Topology,idx);
    %         plot(H,'EdgeLabel',H.Edges.CQI);

            Graph_labels = {['EndNodes1'], ['EndNodes2'], ['CQI'], ['Capacity']};
            % IAB data
            EndNodes1 = H.Edges.EndNodes(:,1);
            EndNodes2 = H.Edges.EndNodes(:,2);
            Graph_data_matrix = zeros(length(EndNodes1)+1,4);
            Graph_data_matrix(:,1) = [EndNodes1; 0];
            Graph_data_matrix(:,2) = [EndNodes2; 0];
            Graph_data_matrix(:,3) = [H.Edges.CQI; 0];
            Graph_data_matrix(:,4) = [H.Edges.Capacity; 0];

            % Donor data
            Graph_data_matrix(end,1) = sum(net_obj.Topology.Edges.Capacity(idxDL));
            Graph_data_matrix(end,2) = floor(mean(net_obj.Topology.Edges.CQI(idxDL)));
            Graph_data_matrix(end,3) = sum(net_obj.Topology.Edges.Capacity(idxUL));
            Graph_data_matrix(end,4) = floor(mean(net_obj.Topology.Edges.CQI(idxUL)));


            % Convert database matrix to table
            Graph_database_Table = table(Graph_data_matrix);
            Graph_database_Table = splitvars(Graph_database_Table);
            Graph_database_Table.Properties.VariableNames = Graph_labels;

            % Save graph data
            filename = ['iteration', num2str(iter)];
            writetable(Graph_database_Table,['IAB-GraphData\',filename,'.csv'])
        end
    end
end