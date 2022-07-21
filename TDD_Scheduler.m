function [TDD_stat] = TDD_Scheduler(timestep,method,disp_flag)
% TDD_Scheduler get current time-step and sheduling methos for TDD and
% return th current state: Uplink/Downlink.
% set disp_flag to 1 for print TDD state for current time-step. 
% set disp_flag to 0 for diseble printing.  
    switch method
        
        case 'scheduler_01'
            schduling = ['U' 'D'];
            index = mod(timestep, length(schduling));
            index(index==0) = length(schduling);
            TDD_stat = schduling(index);
            
            if TDD_stat == 'U' && disp_flag
                disp(['timestep = ',num2str(timestep),' : Uplink Slot'])
            elseif TDD_stat == 'D' && disp_flag
                disp(['timestep = ',num2str(timestep),' : Downlink Slot'])
            end
            
        case 'scheduler_02'
            disp('non')
            TDD_stat = NaN;
    end
end

