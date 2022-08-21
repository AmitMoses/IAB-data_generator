classdef network
   properties
       users
       IABnodes
       IABdonors
       disconnected_users
       Topology
   end
   methods
       %% get_locaions
       function [ue_x_pos, ue_y_pos, BS_x_pos, BS_y_pos, IAB_x_pos, IAB_y_pos] = get_locaions(obj)
           [ue_x_pos] = [obj.users(:).x_pos];
           [ue_y_pos] = [obj.users(:).y_pos];
           [BS_x_pos] = [obj.IABdonors(:).x_pos];
           [BS_y_pos] = [obj.IABdonors(:).y_pos];
           [IAB_x_pos] = [obj.IABnodes(:).x_pos];
           [IAB_y_pos] = [obj.IABnodes(:).y_pos];
       end
       
       %% Calc_distance
       function Distance_Matrix = Calc_distance(obj,x1,y1,x2,y2)
            L1 = length(x1);
            L2 = length(x2);
            BS_x_pos_span = repmat(x2,1,L1)';
            ue_x_pos_span = repmat(x1,L2,1); ue_x_pos_span = ue_x_pos_span(:);
            BS_y_pos_span = repmat(y2,1,L1)';
            ue_y_pos_span = repmat(y1,L2,1); ue_y_pos_span = ue_y_pos_span(:);
            % distance = sqrt((x1-x2)^2 + (y1-y2)^2) @ 2 point (x1,y1),(x2,y2)
            Distance_Vector = sqrt((ue_x_pos_span - BS_x_pos_span).^2 + (ue_y_pos_span - BS_y_pos_span).^2);
            Distance_Matrix = reshape(Distance_Vector,L2,L1);     % col=(x1,y1) row=(x2,y2)
            Distance_Matrix(Distance_Matrix==0) = NaN;
       end
       
       %% RFloss
       function L = RFloss(obj,Distance_km, freq_MHz,method)
            switch method
                case 'Free-Space'
                    L = 32.44 + 20*log10(Distance_km) + 20*log10(freq_MHz);
            end     
       end
       
       %% Pathloss_calculation
       function [PathLoss_Matrix_ue2allBS] = Pathloss_calculation(obj,method)
           [ ue_x_pos, ue_y_pos, BS_x_pos, BS_y_pos, IAB_x_pos, IAB_y_pos] = get_locaions(obj);
           Distance_Matrix_ue2BS = Calc_distance(obj,...
               [ue_x_pos IAB_x_pos],[ue_y_pos IAB_y_pos],...
               [IAB_x_pos BS_x_pos],[IAB_y_pos BS_y_pos]);
           Distance_Matrix_Km_ue2BS = Distance_Matrix_ue2BS/1e3;
           frequncy_matrix_MHz_ue2BS = repmat([obj.IABdonors(:).freq obj.IABnodes(:).freq].'./1e6 , 1, length(obj.users) + length(obj.IABnodes) );
           PathLoss_Matrix_ue2allBS = RFloss(obj,Distance_Matrix_Km_ue2BS, frequncy_matrix_MHz_ue2BS,method);
%            switch method
%                case 'Free-Space'
%                    PathLoss_Matrix_ue2allBS = 32.44 + 20*log10(Distance_Matrix_Km_ue2BS) + 20*log10(frequncy_matrix_MHz_ue2BS);
%                case 'ABG'
%                    alpha = 3.3;
%                    beta = 17.6;
%                    gamma = 2;
%                    PathLoss_Matrix_ue2allBS = 10*alpha*log10(Distance_Matrix_ue2BS) + beta + 10*gamma*log10(frequncy_matrix_MHz_ue2BS./1e3);
%            end
       end
       
       %% plot_network_location
       function plot_network_location(obj)
           disp('plot all units location');
           [ ue_x_pos, ue_y_pos, BS_x_pos, BS_y_pos, IAB_x_pos, IAB_y_pos] = get_locaions(obj);
           figure()
           scatter(ue_x_pos, ue_y_pos,'*') % print UE locations
           hold on
           scatter(BS_x_pos, BS_y_pos,'s') % print gNB locations
           hold on
           scatter(IAB_x_pos, IAB_y_pos,'o') % print IAB locations
           grid on
           grid minor
           legend('UE','IAB donor','IAB Node')
           xlabel('[meters]')
           ylabel('[meters]')
           
           % print UE ID on plot
           hold on
           txt = num2cell([obj.users(:).ID]);
           dx = 5;
           text(ue_x_pos + dx, ue_y_pos,txt)
           
           % print BS ID on plot
           hold on
           txt = num2cell([obj.IABdonors(:).ID]);
           text(BS_x_pos + dx, BS_y_pos,txt)
           
           % print IAB ID on plot
           hold on
           txt = num2cell([obj.IABnodes(:).ID]);
           text(IAB_x_pos + dx, IAB_y_pos,txt)
       end
       
       %% Generate Graph Object representation to the network topology 
       function [G, G_id] = network_topology(obj)
           % the the network and generate Graph Object represention to the
           % network.
           % Output:    G - Graph Object with units names for plots
           %            G_id - Graph Object with information for
           %            calculation
           global CQI2SNR
           
           % Generate Graph Object
           % s, t & w are Graph varibles
           s = [];
           t = [];
           Weight = [];
           
           % Downlink connections 
           index = 1;
           while index <= length(obj.IABdonors)
               connected_UE = [obj.IABdonors(index).UE_con];
               Service_BS = repelem([obj.IABdonors(index).ID], length(connected_UE));
               CQI_DL = obj.IABdonors(index).UE_con_CQI;
               s = [s connected_UE];
               t = [t Service_BS];
               Weight = [Weight CQI_DL];
               index = index +1;
           end
           
           % IAB DL + UL connections 
           index = 1;
           while index <= length(obj.IABnodes)
               
               connected_UE = [obj.IABnodes(index).gNB.UE_con];
               Service_BS = repelem([obj.IABnodes(index).gNB.ID], length(connected_UE));
               
               connected_BS = [obj.IABnodes(index).UE.BS_con_id];
               Service_UE = [obj.IABnodes(index).UE.ID];
               Service_UE = repelem(Service_UE,length(connected_BS));
               
               CQI_DL = obj.IABnodes(index).gNB.UE_con_CQI;
               CQI_UL = obj.IABnodes(index).UE.BS_con_CQI;
               s = [s connected_UE connected_BS];
               t = [t Service_BS Service_UE];
               Weight = [Weight CQI_DL CQI_UL];
               index = index +1;
           end
           
           % Uplink connections
           index = 1; 
           Service_BS = [];
           while index <= length(obj.users)
               Service_BS = [obj.users(index).BS_con_id];
               CQI_DL = obj.users(index).BS_con_CQI;
               s = [s Service_BS];
               t = [t obj.users(index).ID];
               Weight = [Weight obj.users(index).BS_con_CQI];
               index = index +1;
           end
           
           if length(s) > length(Weight)
               return 
           end
           
           % Change id number to String name
           [~ , s_idx_UE] = find(s <= length(obj.users));
           s_val_UE = s(s_idx_UE);
           s_val_UE_str = [repmat('UE_{',length(s_val_UE),1), num2str(s_val_UE'),repmat('}',length(s_val_UE),1)];
           s_val_UE_str = cellstr(s_val_UE_str);
           
           [~ , s_idx_BS] = find(s > length(obj.users) + length(obj.IABnodes));
           s_val_BS = s(s_idx_BS);
           s_val_BS_str = [repmat('BS_{',length(s_val_BS),1), num2str(s_val_BS'),repmat('}',length(s_val_BS),1)];
           s_val_BS_str = cellstr(s_val_BS_str);
           
           [~ , s_idx_IAB] = find(s > length(obj.users) & s <= length(obj.IABnodes)+length(obj.users));
           s_val_IAB = s(s_idx_IAB);
           s_val_IAB_str = [repmat('IAB_{',length(s_val_IAB),1), num2str(s_val_IAB'),repmat('}',length(s_val_IAB),1)];
           s_val_IAB_str = cellstr(s_val_IAB_str);
           
           s_str = cell(size(s));
           s_str(s_idx_UE) = s_val_UE_str;
           s_str(s_idx_BS) = s_val_BS_str;
           s_str(s_idx_IAB) = s_val_IAB_str;
           
           
           [~ , t_idx_UE] = find(t <= length(obj.users));
           t_val_UE = t(t_idx_UE);
           t_val_UE_str = [repmat('UE_{',length(t_val_UE),1), num2str(t_val_UE'),repmat('}',length(t_val_UE),1)];
           t_val_UE_str = cellstr(t_val_UE_str);
           
           [~ , t_idx_BS] = find(t > length(obj.users) + length(obj.IABnodes));
           t_val_BS = t(t_idx_BS);
           t_val_BS_str = [repmat('BS_{',length(t_val_BS),1), num2str(t_val_BS'),repmat('}',length(t_val_BS),1)];
           t_val_BS_str = cellstr(t_val_BS_str);
           
           [~ , t_idx_IAB] = find(t > length(obj.users) & t <= length(obj.IABnodes) + length(obj.users));
           t_val_IAB = t(t_idx_IAB);
           t_val_IAB_str = [repmat('IAB_{',length(t_val_IAB),1), num2str(t_val_IAB'),repmat('}',length(t_val_IAB),1)];
           t_val_IAB_str = cellstr(t_val_IAB_str);
           
           t_str = cell(size(t));
           t_str(t_idx_UE) = t_val_UE_str;
           t_str(t_idx_BS) = t_val_BS_str;
           t_str(t_idx_IAB) = t_val_IAB_str;
           
           % Gnerate graph
           G = digraph(s_str,t_str,Weight);
           G_id = digraph(s,t,Weight);
           G_id.Edges.CQI = G_id.Edges.Weight;
           G_id.Edges.SNR_delay = CQI2SNR(end,1) - G_id.Edges.Weight + 1;
           G_id.Edges.queue_delay = zeros(length(G_id.Edges.EndNodes),1);
%            net.Topology.Edges.SNR_delay = CQI2SNR(net.Topology.Edges.CQI+1,2);
%            G_id.Edges.Weight = CQI2SNR(end,1) - G_id.Edges.Weight + 1;
           G_id.Edges.Weight = G_id.Edges.SNR_delay + G_id.Edges.queue_delay;
           G_id.Edges.Capacity = zeros(size(G_id.Edges.Weight));
           
           % locations node features
           [ ue_x_pos, ue_y_pos, BS_x_pos, BS_y_pos, IAB_x_pos, IAB_y_pos] = get_locaions(obj);
           Xloc = [ue_x_pos, IAB_x_pos, BS_x_pos];
           Yloc = [ue_y_pos, IAB_y_pos, BS_y_pos];
           G.Nodes.Xloc = Xloc.';
           G.Nodes.Yloc = Yloc.';
           
       end
       
       %% Plot network topology
       function obj = plot_network_topology(obj)
           disp('Plot Network Topology');
           [G,~] = network_topology(obj);
           figure()
%            plot(G,'EdgeLabel',G.Edges.Weight, 'XData',G.Nodes.Xloc, 'YData', G.Nodes.Yloc);
           plot(G,'EdgeLabel',G.Edges.Weight);
           
       end
             
       %% Data generation
       function DataMatrix = DataGenerator(obj,UnitNum)
           % Generate DataMatrix which contaion the sore and destanation of
           % each data in the format [Senders,Recivers]
           % Senders = [Ue_Num IABnode_num IABdonor_Num]
           % Recivers = [Ue_Num IABnode_num IABdonor_Num]
           
            global Ue_Num
            global IABnode_num
            global IABdonor_Num
            DataMatrix = zeros(UnitNum);
            % Uplink
            DataMatrix(floor(1:Ue_Num/4),end) = 5;
            DataMatrix(floor(Ue_Num/4+1:Ue_Num),end) = 2.5;
            % Downlink
            DataMatrix(end,floor(1:Ue_Num/4)) = 2;
            DataMatrix(end,floor(Ue_Num/4+1:Ue_Num)) = 0.5;
       end
       
       
       %% branch capacity calculation
       function obj = Datapath(obj,UnitNum)
           % This function find all the paths from the senders to recivers
           % and update the load on the Graph Object on the Topology field 
           % of the network. The data that each user send is equal to all
           % recivers, so the function use 'Minimum Spanning Tree' for load
           % calculation

            global IABnode_num
            global IABdonor_Num
            DataMatrix = DataGenerator(obj,UnitNum);
            [Senders, Recivers] = find(DataMatrix);
            data_capacity_matrix = zeros(size(obj.Topology.Edges.EndNodes,1),UnitNum);
            allDataPaths  = NaN(length(Senders),IABdonor_Num + IABnode_num);
            
            % Delete overflow nodes -- old version note that its need improvment 
            Delete_Senders = find(Senders>size(obj.Topology.Nodes,1));
            Delete_Recivers = find(Recivers>size(obj.Topology.Nodes,1));
            Delete_vector = ([Delete_Senders ;Delete_Recivers]);
            Senders(Delete_vector) = [];
            Recivers(Delete_vector) = [];
            
            for link = 1:length(Senders)
                % shotrstpath:      'unweighted' - Breadth-First compution
                %                   that treat all edge weights as 1
                %                   'positive' - use Dijkstra algorithem to
                %                   find the shortest path
               path = shortestpath(obj.Topology, Senders(link), Recivers(link),'Method','positive');
               allDataPaths(link, 1:length(path)) = path;
               
               % Get a path and place the load (data) on all the roads in
               % the Graph Object. Make a node-to-node matrix that each row
               % represent a raod in the path 
               path(isempty(path)) = 0; % if there is no connection set 0 (insted of [])
               x1 = repelem(path,2); % replicate the vector elements
               if length(x1(2:end-1)) > 2
                    v2 = reshape(x1(2:end-1),2,[]).';
               else
                    v2 = x1(2:end-1);
               end
%                sort_edges = sort(reshape(v2,[],2),1);
               sort_edges = v2;
               ai = find(ismember(obj.Topology.Edges.EndNodes, sort_edges,'rows'));
%                [~,ai,~] = intersect(obj.Topology.Edges.EndNodes,sort_edges,'rows');
%                data_capacity_matrix(ai, Senders(link)) = data_capacity_matrix(ai, Senders(link)) + DataMatrix(Senders(link),Recivers(link));
               obj.Topology.Edges.Capacity(ai) = obj.Topology.Edges.Capacity(ai) + DataMatrix(Senders(link),Recivers(link));
            end
%             obj.Topology.Edges.Capacity = obj.Topology.Edges.Capacity + Load*sum(double(data_capacity_matrix > 0),2); 
%             obj.Topology.Edges.Capacity = obj.Topology.Edges.Capacity + sum(double(data_capacity_matrix),2); 
       end
       
       %% Update wights
       function obj = update_weights(obj)
           obj.Topology.Edges.Weight = obj.Topology.Edges.SNR_delay + obj.Topology.Edges.queue_delay;
       end
       
       %% branch capacity calculation
       function obj = Random_Datapath(obj,UnitNum)
           % This function find all the paths from the senders to recivers
           % and update the load on the Graph Object on the Topology field 
           % of the network. The data that each user send is equal to all
           % recivers, so the function use 'Minimum Spanning Tree' for load
           % calculation

            global IABnode_num
            global IABdonor_Num
            DataMatrix = DataGenerator(obj,UnitNum);
            [Senders, Recivers] = find(DataMatrix);
            p = randperm(length(Senders));  % permutation vector
            Senders = Senders(p);
            Recivers = Recivers(p);
%             data_capacity_matrix = zeros(size(obj.Topology.Edges.EndNodes,1),UnitNum);
            allDataPaths  = NaN(length(Senders),IABdonor_Num + IABnode_num);
            
            % Delete overflow nodes -- old version note that its need improvment 
            Delete_Senders = find(Senders>size(obj.Topology.Nodes,1));
            Delete_Recivers = find(Recivers>size(obj.Topology.Nodes,1));
            Delete_vector = ([Delete_Senders ;Delete_Recivers]);
            Senders(Delete_vector) = [];
            Recivers(Delete_vector) = [];
            
            for link = 1:length(Senders)
                % shotrstpath:      'unweighted' - Breadth-First compution
                %                   that treat all edge weights as 1
                %                   'positive' - use Dijkstra algorithem to
                %                   find the shortest path
               path = shortestpath(obj.Topology, Senders(link), Recivers(link),'Method','positive');
               allDataPaths(link, 1:length(path)) = path;
               
               % Get a path and place the load (data) on all the roads in
               % the Graph Object. Make a node-to-node matrix that each row
               % represent a raod in the path 
               path(isempty(path)) = 0; % if there is no connection set 0 (insted of [])
               x1 = repelem(path,2); % replicate the vector elements
               if length(x1(2:end-1)) > 2
                    v2 = reshape(x1(2:end-1),2,[]).';
               else
                    v2 = x1(2:end-1);
               end
               sort_edges = v2;
               ai = find(ismember(obj.Topology.Edges.EndNodes, sort_edges,'rows'));
               obj.Topology.Edges.Capacity(ai) = obj.Topology.Edges.Capacity(ai) + DataMatrix(Senders(link),Recivers(link));
               
               % add queue
               allocate_BW = obj.Topology.Edges.allocate_BW(ai);  %Hz TotalBW=550e6, num_of_links=110*2
               RB_BW = 180e3;       % Hz
               num_RB = allocate_BW/RB_BW;
               packet_size = 100e3; %[bps]
               massage_size = obj.Topology.Edges.Capacity(ai) * 1e6; % [bps] 
               queue_length = massage_size / packet_size;
               obj.Topology.Edges.queue_delay(ai) = obj.Topology.Edges.queue_delay(ai) + queue_length./num_RB;
               obj = update_weights(obj);
            end
       end
       
       %% set resource allocation
       function obj = set_resource_allocation(obj, method)
           switch method
               case 'fair'
                   allocate_BW = 2.5e6;  %Hz TotalBW=550e6, num_of_links=110*2
                   allocate_BW = repmat(allocate_BW, size(obj.Topology.Edges(:,1),1),1);
               case 'load'
                   disp('load')
           end
           obj.Topology.Edges.allocate_BW = allocate_BW;
       end
   end
end







