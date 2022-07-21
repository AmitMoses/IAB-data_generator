classdef WaveForm
   properties
       % Bandwidth of the WaveForm
       % FR1 support maximum bandwidth of 100MHZ:
       %    BW_FR1 = [ 5; 10 ; 15; 20; 25; 40; 50; 60; 80; 100 ] MHz;
       % FR2 support maximum bandwidth of 400MHZ:
       %    BW_FR2 = [ 50; 100; 200; 400 ] MHz;
       BW                           % [Hz]
       
       % meu take value from 0 to 5 accurding to FR1/FR2
       % for FR1: meu = 0/1
       % for FR2: meu = 2/3/4 and 5 for Real.16
       meu
       
       % Sub Carrier Spacing is Injective function of meu:
       %    SubCarrierSpace = 2^meu*15KHz
       % Will set up in the WaveForm SetUp Function.
       SubCarrierSpace              % [Hz]
       
       % Duration Time of one Frame.
       % Usally 10 times of one Sub-Frame in secounds.
       % Will set up in the WaveForm SetUp Function.
       Time_Frame                   % [sec]     
       
       % Duration Time of one Sub-Fram in secounds.
       % A default Value
       Time_SubFrame        = 1e-3; % [sec]      
       
       % Duratin Time of one slot in secounds.
       % Will set up in the WaveForm SetUp Function accurding Time_SubFrame
       % and meu
       Time_Slot                    % [sec]
       
       % Number of Slot per one Sub-Frame.
       % The divertion between Time_SubFrame and Time_Slot.
       % It is also represent the Physical Reasoure Block Duration.
       NumSlotPerSubframe
       
       % Number of OFDM symbole per one Time-Slot.
       % A default Value.
       OFDMsymPerSlot       = 14;
       
       % Total Number of Physical Reasoure Block.
       NumPRB
       
       % Number of Subcarriers in one Physical Reasoure Block.
       % A default Value.
       PRB_subcarrier       = 12;
       
   end
   
   methods
      
       function obj = setup(obj, BW, meu)
           obj.BW = BW;
           obj.meu = meu;
           obj.SubCarrierSpace = 2.^meu*15e3;
           obj.Time_Frame = 10*obj.Time_SubFrame;
           obj.Time_Slot = obj.Time_SubFrame / (2^meu);
           obj.NumSlotPerSubframe = obj.Time_SubFrame / obj.Time_Slot;
           obj.NumPRB = floor( BW / (obj.PRB_subcarrier*obj.SubCarrierSpace) );
       end
       
   end
end