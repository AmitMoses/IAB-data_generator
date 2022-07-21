close all
clear 
clc

% FR1 support meu = 0/1 OR Subcarrier spacing of 15/30 KHz
% FR2 support meu = 2/3/4 OR Subcarrier spacing of 60/120/240 KHz
meu = (0:5)';
SubSpace = 2.^meu*15e3;
% FR1 support maximum bandwidth of 100MHZ
% FR2 support maximum bandwidth of 400MHZ
BW_FR1 = [ 5; 10 ; 15; 20; 25; 40; 50; 60; 80; 100]*1e6;
BW_FR2 = [50; 100; 200; 400]*1e6;

% Frame Structure
TimeSubframe = 1e-3;
TimeFrame = 10*TimeSubframe;
OFDMsymPerSlot = 14;    % Number of OFDM symbols per time slot;
slot_lenght = TimeSubframe/(2^meu(1)); 
NumSlotPerSubframe =  TimeSubframe / slot_lenght;

%Bandwidth calculation
N_RB = 28;
PRB_subcarrier = 12;
BW_1RB = PRB_subcarrier * SubSpace(1);
BW_tot = BW_1RB*N_RB;
% input(BW,Sun Carrier Spacing OR meu) --> Output(#PRB,slot_lenght,NumSlotPerSubframe)
