%%*********************************************************
%% Demo file with minimal example for SBSOS
%%
%% Author: Tillmann Weisser
%%*********************************************************

clear all;
clc;

% pop.F = [2 0 0 1; 1 1 0 -1; 0 2 0 1; 0 1 1 -1; 0 0 2 1];
% pop.n = 3;
% pop.I = {1:2,2:3};
% 
% pop.G{1} = [1 0 0 1];
% pop.G{2} = [0 1 0 1];
% pop.G{3} = [0 0 1 1];
% pop.G{4} = [2 0 0 -1; 0 2 0 -1; 0 0 0 1];
% pop.G{5} = [0 0 2 -1; 0 2 0 -1; 0 0 0 1];
% 
% pop.J = {[1,2,4],[2,3,5]};
pop.F = [2 0  0 1; 1 1 0 -1; 0 2 0 1];
pop.n = 3;
pop.I = {1:3};

pop.G{1} = [1 0 0 -1;0 0 0 10];
pop.G{2} = [0 1 0 -1;0 0 0 5];
pop.G{3} = [0 0 1 -1; 0 0 0 2];
%pop.G{3} = [2 0 0 1; 0 2 0 1; 0 0 0 -1];

pop.J = {[3]};

pop.k=1; 
pop.d=1;

sdp = gendata2(pop,'SBSOS');

sol = csol(sdp,'sdpt3');

psol = postproc(pop,sdp,sol);

