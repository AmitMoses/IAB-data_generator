classdef saving
    properties
        ue_data
        iab_data
    end
    methods
        %%
        function obj = set(obj, UE_database, IAB_database)
            obj.ue_data = UE_database;
            obj.iab_data = IAB_database;
        end
        
        %%
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
        
        %%
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
    end
end