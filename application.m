classdef application
    properties
        onoff   % on - 1 / off - 0
        rate    % Data rate transsmit [Mbps]
        dest    % IP of destanation
    end
    methods
        function obj = off(obj)
            obj.onoff = 0;
            obj.rate = [];
            obj.dest = [];
        end
        function obj = on(obj, rate, dest)
            obj.onoff = 1;
            obj.rate = rate;
            obj.dest = dest;
        end
    end
end