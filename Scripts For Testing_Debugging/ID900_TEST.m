% ID900 Connect Test
clear;
clc;
% Object with all required functions
MyFuncs = functionsContainer;
MyFuncs.AddPathFunc("LAB"); 
tc = py.ID900_Func.init_ID900();

photon_count = py.ID900_Func.query_photon_counter(tc); 
