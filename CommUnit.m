classdef CommUnit
    properties
        ID                  % ID name
        x_pos               % x axis position
        y_pos               % y axis position
        freq                % Tx/Rx frequncy
        Ptx                 % Transmit power [dBm]
        app = application
    end
    
    methods
        
        function obj = set_app(obj, onoff, rate, dest)
            if onoff == 1
                obj.app = on(obj.app, rate, dest);
            elseif onoff == 0
                obj.app = off(obj.app);
            end
        end
        
    end
end