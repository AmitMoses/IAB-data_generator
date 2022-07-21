close all
clear
clc

global Ue_Num
global IABnode_num
global IABdonor_Num
global CQI2SNR

%% Flags
hop_1_IAB = 0;
GenerateDatabase = 0; 
GenerateGraphData = 0;

%% Iteration's number
% MaxIterations = 10000;
% MaxIterations = 1000;
MaxIterations = 1;

%% Simulation Parameters:
IABnode_num = 9;
IABdonor_Num = 1;
Ue_Num = (IABnode_num+IABdonor_Num)*10;
% Ue_Num = (IABnode_num+IABdonor_Num);
UnitNum = Ue_Num + IABnode_num + IABdonor_Num;
BS_frequncy = [ 3.9e9 ]*( ones(1,IABdonor_Num+IABnode_num) ) ;
Total_Bandwith = 10e6; %[MHz]
AreaSize = (40000);
% data : AreaSize = (10000);
% data_V2 : AreaSize = (20000);
% data_V3 : AreaSize = (40000);

CQI2SNR = [
            0 -inf
            1 -6.82;
            2 -3.44;
            3 0.53;
            4 3.79;
            5 5.8;
            6 8.08;
            7 9.76;
            8 11.72;
            9 13.49;
            10 15.87;
            11 17.73;
            12 19.50;
            13 21.32;
            14 23.51;
            15 25;
            15 inf];



%% Network Calibration
UE_database = zeros(MaxIterations,(IABnode_num+IABdonor_Num)*5);
IAB_database = zeros(MaxIterations,(IABnode_num+IABdonor_Num)*5);
tic
for iter = 1:MaxIterations

    net = network;
    net.users = UE;
    net.IABnodes = IABnode;
    net.IABdonors = IABdonor;

    % UE calibration
    for UE_idx=1:Ue_Num
        net.users(UE_idx).ID = UE_idx;
        net.users(UE_idx).x_pos = randi(AreaSize);   
        net.users(UE_idx).y_pos = randi(AreaSize);
        net.users(UE_idx).Ptx = 23;
        net.users(UE_idx).BW = Total_Bandwith;
    end

    % IAB calibration
    for IAB_idx = 1:IABnode_num
        net.IABnodes(IAB_idx) = IABnode;
        net.IABnodes(IAB_idx) = set_IABnode(net.IABnodes(IAB_idx),...
            IAB_idx + Ue_Num,...
            randi(AreaSize),...
            randi(AreaSize),...
            BS_frequncy(IAB_idx + IABdonor_Num),...
            30,...
            Total_Bandwith,...
            Total_Bandwith);
    end

    % BS calibration
    for BS_idx=1:IABdonor_Num
        net.IABdonors(BS_idx).ID = BS_idx + Ue_Num + IABnode_num;      
        net.IABdonors(BS_idx).x_pos  = randi(AreaSize); 
        net.IABdonors(BS_idx).y_pos = randi(AreaSize);   
        net.IABdonors(BS_idx).freq = BS_frequncy(BS_idx);
        net.IABdonors(BS_idx).Ptx = 30;
        net.IABdonors(BS_idx).BW = Total_Bandwith;
    end

    

    %% Network connections
    % calculate Path Loss between all users to BS (Network object method)
    % [Pathloss_Matrix_ue2allBS] = Pathloss_calculation(net,'ABG');
    [Pathloss_Matrix_ue2allBS] = Pathloss_calculation(net,'Free-Space');
    Pathloss_Matrix_ue2BS = Pathloss_Matrix_ue2allBS(IABnode_num+IABdonor_Num:end,1:Ue_Num);
    Pathloss_Matrix_ue2IAB = Pathloss_Matrix_ue2allBS(1:IABnode_num,1:Ue_Num);
    Pathloss_Matrix_IAB2BS = Pathloss_Matrix_ue2allBS(IABnode_num+IABdonor_Num:end,Ue_Num+1:end);
    Pathloss_Matrix_IAB2IAB = Pathloss_Matrix_ue2allBS(1:IABnode_num , Ue_Num+1:end);

    

    % IAB search for best BS to connect
    for UE_idx=1:IABnode_num
        % IAB establish connection with BS
        if hop_1_IAB == 1
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, Pathloss_Matrix_IAB2BS(:,UE_idx), inf ,1); 
        else
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, Pathloss_Matrix_IAB2BS(:,UE_idx), Pathloss_Matrix_IAB2IAB(:,UE_idx) ,1); 
        end
    end

    % UE search for best BS to connect
    for UE_idx=1:Ue_Num
        % UE establish connection with BS
        [net] =...
            Connect2BS(net.users(UE_idx), net, Pathloss_Matrix_ue2BS(:,UE_idx), Pathloss_Matrix_ue2IAB(:,UE_idx),0); 
    end
    
    % BS establish connection with UE
    for BS_idx1 = 1:IABdonor_Num
            net.IABdonors(BS_idx1) =...
            Connect2UE(net.IABdonors(BS_idx1),[net.users net.IABnodes.UE] , Pathloss_Matrix_ue2allBS(BS_idx1 + length(net.IABnodes),:),CQI2SNR);
    end

    % IAB establish connection with UE
    for IAB_idx = 1:IABnode_num
        net.IABnodes(IAB_idx).gNB =...
            Connect2UE(net.IABnodes(IAB_idx).gNB, [net.users net.IABnodes.UE] , Pathloss_Matrix_ue2allBS(IAB_idx,:),CQI2SNR);
    end
    
    
    % Plot all units locations & Network Topology
    if MaxIterations == 1
        plot_network_location(net)
        plot_network_topology(net);
    end
    
    %% Data Transfer
    [~, net.Topology] = network_topology(net);
    for UE_idx=1:Ue_Num
        % UE hops from BS
        path_to_BS = shortestpath(net.Topology, net.users(UE_idx).ID ,UnitNum,'Method','positive');
        hops = length(path_to_BS) - 1;
        net.users(UE_idx).hops = hops;
    end
    net = Random_Datapath(net, UnitNum);
%     net = Datapath(net, UnitNum);
    
    


    % %% TDD schduling
    % disp_flag = 1;
    % for timestep = 1:10
    %     [TDD_slot] = TDD_Scheduler(timestep,'scheduler_01',disp_flag);
    % end
    % 
    % resource_allocation(net.IABdonors(1));

    %% Save data
    link_size = size(net.Topology.Edges.EndNodes,1);
    % saving UE database into UE_database matrix
    for i=0:Ue_Num-1
        idxDL = find(net.Topology.Edges.EndNodes(:,1) == i+1);
        idxUL = find(net.Topology.Edges.EndNodes(:,2) == i+1);
        UE_database(iter,i*5+1) = net.users(i+1).BS_con_id;               % connected BS
        UE_database(iter,i*5+2) = net.Topology.Edges.Capacity(idxDL);     % DL_app
        UE_database(iter,i*5+3) = net.Topology.Edges.CQI(idxDL);          % DL_CQI
        UE_database(iter,i*5+4) = net.Topology.Edges.Capacity(idxUL);     % UL_app
        UE_database(iter,i*5+5) = net.Topology.Edges.CQI(idxUL);          % UL_CQI
    end
    
    % saving IABnodes database into IAB_database matrix
    for i=0:IABnode_num -1
         IABnode_ID = net.IABnodes(i+1).ID;
         connected_gNB = net.IABnodes(i+1).UE.BS_con_id;

%          idxDL = intersect(...
%             find(net.Topology.Edges.EndNodes(:,1) > Ue_Num ),...
%             find(net.Topology.Edges.EndNodes(:,2) == net.IABnodes(i+1).ID));
%          
%          idxUL = intersect(...
%             find(net.Topology.Edges.EndNodes(:,2) > Ue_Num ),...
%             find(net.Topology.Edges.EndNodes(:,1) == net.IABnodes(i+1).ID));
        
         idxDL = find(ismember(net.Topology.Edges.EndNodes, [connected_gNB, IABnode_ID],'rows'));
         idxUL = find(ismember(net.Topology.Edges.EndNodes, [IABnode_ID, connected_gNB],'rows'));
         
         
        if isnan(net.IABnodes(i+1).UE.BS_con_id) 
            IAB_database(iter,i*5+1) = 0;
        else
            IAB_database(iter,i*5+1) = net.IABnodes(i+1).UE.BS_con_id;          % Connected BS
        end
        
        if isnan(net.Topology.Edges.Capacity(idxDL))
            IAB_database(iter,i*5+2) = 0;
        else
            IAB_database(iter,i*5+2) = net.Topology.Edges.Capacity(idxDL);      % DL_app
        end
        
        if isnan(net.Topology.Edges.CQI(idxDL))
            IAB_database(iter,i*5+3) = 0;
        else
            IAB_database(iter,i*5+3) = net.Topology.Edges.CQI(idxDL);           % DL_CQI
        end
        
        if isnan(net.Topology.Edges.Capacity(idxUL))
            IAB_database(iter,i*5+4) = 0;
        else
            IAB_database(iter,i*5+4) = net.Topology.Edges.Capacity(idxUL);      % UL_app
        end
        
        if isnan(net.Topology.Edges.CQI(idxUL))
            IAB_database(iter,i*5+5) = 0;
        else
            IAB_database(iter,i*5+5) = net.Topology.Edges.CQI(idxUL);           % UL_CQI
        end
    end
    
    % saving IABdonors database into IAB_database matrix
    for i=0:IABdonor_Num -1
         idxDL = intersect(...
            find(net.Topology.Edges.EndNodes(:,2) <= Ue_Num ),...
            find(net.Topology.Edges.EndNodes(:,1) == net.IABdonors(i+1).ID));
         idxUL = intersect(...
            find(net.Topology.Edges.EndNodes(:,1) <= Ue_Num ),...
            find(net.Topology.Edges.EndNodes(:,2) == net.IABdonors(i+1).ID));
        IAB_database(iter,(i + IABnode_num)*5+1) = -1;          % Connected BS
        if isnan(sum(net.Topology.Edges.Capacity(idxDL))) 
            IAB_database(iter,(i + IABnode_num)*5+2) = 0;
        else
            IAB_database(iter,(i + IABnode_num)*5+2) = sum(net.Topology.Edges.Capacity(idxDL));         % DL_app
        end
        
        if isnan(floor(mean(net.Topology.Edges.CQI(idxDL))))
            IAB_database(iter,(i + IABnode_num)*5+3) = 0;
        else
            IAB_database(iter,(i + IABnode_num)*5+3) = floor(mean(net.Topology.Edges.CQI(idxDL)));      % DL_CQI
        end
        
        if isnan(sum(net.Topology.Edges.Capacity(idxUL)))
            IAB_database(iter,(i + IABnode_num)*5+4) = 0;
        else
            IAB_database(iter,(i + IABnode_num)*5+4) = sum(net.Topology.Edges.Capacity(idxUL));         % UL_app
        end
        
        if isnan(floor(mean(net.Topology.Edges.CQI(idxUL))))
            IAB_database(iter,(i + IABnode_num)*5+5) = 0;
        else
            IAB_database(iter,(i + IABnode_num)*5+5) = floor(mean(net.Topology.Edges.CQI(idxUL)));      % UL_CQI
        end
    end
    
    if mod(iter,10) == 0
        disp(['complete iteration: ', num2str(iter),', time = ',num2str(toc)] )
    end
    
    % Saving IAB-nodes Graph Data
    if GenerateGraphData
        idx =(1:IABnode_num+IABdonor_Num)+Ue_Num;
        H = subgraph(net.Topology,idx);
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
        Graph_data_matrix(end,1) = sum(net.Topology.Edges.Capacity(idxDL));
        Graph_data_matrix(end,2) = floor(mean(net.Topology.Edges.CQI(idxDL)));
        Graph_data_matrix(end,3) = sum(net.Topology.Edges.Capacity(idxUL));
        Graph_data_matrix(end,4) = floor(mean(net.Topology.Edges.CQI(idxUL)));
        

        % Convert database matrix to table
        Graph_database_Table = table(Graph_data_matrix);
        Graph_database_Table = splitvars(Graph_database_Table);
        Graph_database_Table.Properties.VariableNames = Graph_labels;

        % Save graph data
        filename = ['iteration', num2str(iter)];
        writetable(Graph_database_Table,['IAB-GraphData\',filename,'.csv'])
    end
    
end

%% Save database as table
% Generate UE database labels
if GenerateDatabase
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
    UE_database_Table = table(UE_database);
    UE_database_Table = splitvars(UE_database_Table);
    UE_database_Table.Properties.VariableNames = UE_database_labels;

% Generate IAB database labels
IAB_database_labels = [];
    for  i=1:IABnode_num + IABdonor_Num
        str_0 = ['IAB', num2str(i + Ue_Num),'-Con_BS '];
        str_1 = ['IAB', num2str(i+Ue_Num),'-DL_app '];
        str_2 = ['IAB', num2str(i+Ue_Num),'-DL_CQI '];
        str_3 = ['IAB', num2str(i+Ue_Num),'-UL_app '];
        str_4 = ['IAB', num2str(i+Ue_Num),'-UL_CQI '];
        IAB_database_labels = [IAB_database_labels,str_0, str_1,str_2,str_3,str_4];
    end
    IAB_database_labels = split(IAB_database_labels);
    IAB_database_labels = IAB_database_labels(1:end-1).';

    % Convert database matrix to table
    IAB_database_Table = table(IAB_database);
    IAB_database_Table = splitvars(IAB_database_Table);
    IAB_database_Table.Properties.VariableNames = IAB_database_labels;

    %Save to csv file
    disp('Saving data into csv file...')
    writetable(UE_database_Table,'UE_database.csv')  
    writetable(IAB_database_Table,'IAB_database.csv')  
    disp('Done')

end


sum(net.Topology.Edges.Capacity)
figure()
hist(net.Topology.Edges.CQI)
mean(net.Topology.Edges.CQI)

