classdef IABnode < CommUnit
    properties
        UE = UE
        gNB = IABdonor 
        connect2Donor = 0
    end
    
    methods
       
        %% Sum Apps: 
        %  add all coneccted ue data that connectet throght gNB to UE app
        
        %% Set:
        %  set basic parameters as UE.ID & gNB.ID
        function obj = set_IABnode(obj,ID, x_pos, y_pos, freq, Ptx, ueBW, gnbBW)
            obj.ID = ID;
            obj.x_pos = x_pos;
            obj.y_pos = y_pos;
            obj.freq = freq;
            
            obj.UE.ID = ID;
            obj.UE.x_pos = x_pos;
            obj.UE.y_pos = y_pos;
            obj.UE.freq = freq;
            obj.UE.Ptx = Ptx;
            obj.UE.BW = ueBW;
            
            obj.gNB.ID = ID;
            obj.gNB.x_pos = x_pos;
            obj.gNB.y_pos = y_pos;
            obj.gNB.freq = freq;
            obj.gNB.Ptx = Ptx;
            obj.gNB.BW = gnbBW;
        end
    end
end