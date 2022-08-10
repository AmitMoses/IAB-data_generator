close all
clear
clc

global Ue_Num
global IABnode_num
global IABdonor_Num
global CQI2SNR
global max_bachaul_num
global AreaSize
global Total_Bandwith
global BS_frequncy
global IAB_backhaul_method
global MaxScenarios
global UnitNum


%% Flags
IAB_backhaul_method = 2;    % 0 -   Direct connection to the IAB-Donor
                            % 1 -   Multi hop to the IAB-Donor
                            % 2 -   Mesh conectivity / multipule backhaul
                            %       connections
GenerateDatabase = 1; 
GenerateGraphData = 1;

%% Iteration's number
% MaxScenarios = 10000;
% MaxScenarios = 1000;
MaxScenarios = 10;

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
num_of_features = 5;
max_bachaul_num = 2;
UE_database = zeros(MaxScenarios,(IABnode_num+IABdonor_Num)*num_of_features);
IAB_database = zeros(MaxScenarios,(IABnode_num*max_bachaul_num+IABdonor_Num)*num_of_features);
% Saving Object
save_obj = saving;
% save_obj.ue_data = UE_database;
% save_obj.iab_data = IAB_database;
save_obj = save_obj.set(UE_database, IAB_database);
tic
for scenario = 1:MaxScenarios

    [net] = Create_Random_Network();
    TimeSlotsInScenario = 10;
    
    for slot = 1:TimeSlotsInScenario     
        net_scenario = Random_Datapath(net, UnitNum);
    %     net = Datapath(net, UnitNum);


        %% Save data
        table_row = (scenario-1)*TimeSlotsInScenario + slot;
        % saving UE database into UE_database matrix
        save_obj = save_obj.update_ue(net_scenario, table_row);
        % saving IAB-Nodes database into IAB_database matrix
        save_obj = save_obj.update_iab(net_scenario, table_row);

        % Saving IAB-nodes Graph Data
        if GenerateGraphData
            save_obj.IABgraph(net_scenario, table_row)
        end

        % time display
        if mod(table_row,10) == 0
            disp(['complete iteration: ', num2str(table_row),', time = ',num2str(toc)] )
        end
    end
    
end

%% Save database as table
if GenerateDatabase
    save_obj.UE_Table;
    save_obj.IAB_Table;
end

sum(net.Topology.Edges.Capacity)
figure()
hist(net.Topology.Edges.CQI)
mean(net.Topology.Edges.CQI)

