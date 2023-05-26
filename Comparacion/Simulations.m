% Run this code to create the simulation data for comparision.
% The code ensure the correct naming of the files used to save the
% simulations data. 
% Simulation.slx can be used directly, but care must be put in the name of
% the files in block To File1 and To File2 so it match the selected
% friction model

% Simulation with Dahl friction
set_param('Simulation/To File1','FileName','q1_dahl.mat')
set_param('Simulation/To File2','FileName','q2_dahl.mat')
set_param('Simulation/Constant','Value','1')
sim('Simulation',20)

% Simulation with LuGre friction
set_param('Simulation/To File1','FileName','q1_lugre.mat')
set_param('Simulation/To File2','FileName','q2_lugre.mat')
set_param('Simulation/Constant','Value','2')
sim('Simulation',20)

% Simulation with viscous and Coulomb plus Stribeck effect friction
set_param('Simulation/To File1','FileName','q1_vcs.mat')
set_param('Simulation/To File2','FileName','q2_vcs.mat')
set_param('Simulation/Constant','Value','3')
sim('Simulation',20)

% Simulation w/o frcition
set_param('Simulation/To File1','FileName','q1_wof.mat')
set_param('Simulation/To File2','FileName','q2_wof.mat')
set_param('Simulation/Constant','Value','4')
sim('Simulation',20)