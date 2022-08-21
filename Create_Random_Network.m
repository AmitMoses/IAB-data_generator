function [net] = Create_Random_Network()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    global Ue_Num
    global IABnode_num
    global IABdonor_Num
    global CQI2SNR
    global AreaSize
    global Total_Bandwith
    global BS_frequncy
    global IAB_backhaul_method
    global MaxScenarios
    global UnitNum

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
    [Pathloss_Matrix_ue2iab] = Pathloss_calculation(net,'Free-Space');
    Pathloss_Matrix_ue2donor = Pathloss_Matrix_ue2iab(IABnode_num+IABdonor_Num:end,1:Ue_Num);
    Pathloss_Matrix_ue2node = Pathloss_Matrix_ue2iab(1:IABnode_num,1:Ue_Num);
    Pathloss_Matrix_node2donor = Pathloss_Matrix_ue2iab(IABnode_num+IABdonor_Num:end,Ue_Num+1:end);
    Pathloss_Matrix_node2node = Pathloss_Matrix_ue2iab(1:IABnode_num , Ue_Num+1:end);
    Pathloss_Matrix_iab2iab = Pathloss_Matrix_ue2iab(1:end , Ue_Num+1:end);


    

    % IAB search for best BS to connect
    for UE_idx=1:IABnode_num
        % IAB establish connection with BS
        % UE-class
        if IAB_backhaul_method == 0
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, Pathloss_Matrix_node2donor(:,UE_idx), inf ,1); 
        elseif IAB_backhaul_method == 1
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, Pathloss_Matrix_node2donor(:,UE_idx), Pathloss_Matrix_node2node(:,UE_idx) ,1); 
        elseif IAB_backhaul_method == 2
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, Pathloss_Matrix_node2donor(:,UE_idx), inf ,1); 
%             [net] =...
%             Connect2BS(net.IABnodes(UE_idx).UE, net, inf, Pathloss_Matrix_node2node(:,UE_idx) ,1); 
        end
    end
    
    if IAB_backhaul_method == 2
        for UE_idx=1:IABnode_num
            [net] =...
            Connect2BS(net.IABnodes(UE_idx).UE, net, inf, Pathloss_Matrix_node2node(:,UE_idx) ,1); 
        end
    end


    % UE search for best BS to connect
    for UE_idx=1:Ue_Num
        % UE establish connection with BS
        [net] =...
            Connect2BS(net.users(UE_idx), net, Pathloss_Matrix_ue2donor(:,UE_idx), Pathloss_Matrix_ue2node(:,UE_idx),0); 
    end
    
    % BS establish connection with UE
    for BS_idx1 = 1:IABdonor_Num
            net.IABdonors(BS_idx1) =...
            Connect2UE(net.IABdonors(BS_idx1),[net.users net.IABnodes.UE] , Pathloss_Matrix_ue2iab(BS_idx1 + length(net.IABnodes),:),CQI2SNR);
    end

    % IAB establish connection with UE
    for IAB_idx = 1:IABnode_num
        net.IABnodes(IAB_idx).gNB =...
            Connect2UE(net.IABnodes(IAB_idx).gNB, [net.users net.IABnodes.UE] , Pathloss_Matrix_ue2iab(IAB_idx,:),CQI2SNR);
    end
    
    
    % Plot all units locations & Network Topology
    if MaxScenarios == 1
%         plot_network_location(net)
        plot_network_topology(net);
    end
    
    % Topology
    [~, net.Topology] = network_topology(net);
    for UE_idx=1:Ue_Num
        % UE hops from BS
        path_to_BS = shortestpath(net.Topology, net.users(UE_idx).ID ,UnitNum,'Method','positive');
        hops = length(path_to_BS) - 1;
        net.users(UE_idx).hops = hops;
    end
    
    net = set_resource_allocation(net, 'fair');
end

