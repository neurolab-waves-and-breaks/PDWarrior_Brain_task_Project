function [Const] = storeConst()
 

Const.participants   = {'P01', 'P03','P04','P05','P07','P08','P09','P10','P12','P13','P14','P15','P16','P17','P18','P19','P20','P21'};


Const.exerc = {'preE', 'postE'};
% Const.meds  = {'ON', 'ON-OFF', 'OFF'};
% Const.meds  = {'ON',  'OFF'};
Const.hand  = {'right', 'left'};

Const.cols.lightgreen   = [191, 220, 190] / 255;
Const.cols.midgreen     = [78, 186, 153] / 255;
Const.cols.lightorange  = [247, 211, 163] / 255;
Const.cols.darkorange   = [229, 159, 31] / 255;
Const.cols.lightblue    = [133, 196, 219] / 255;
Const.cols.midblue      = [64, 177, 200] / 255;
Const.cols.darkblue     = [63, 98, 141] / 255;

Const.allCols = [Const.cols.lightorange; Const.cols.darkorange ; ...
                 Const.cols.lightblue; Const.cols.darkblue];

