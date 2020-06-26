diary
diary('mnew90000.txt');


addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/Learn_invariants');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SOSTOOLS-3.300');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/pics');
addpath('/home/research/g.dehghanpoor/libraries/tbxmanager/toolboxes/sedumi/1.3/glnxa64/sedumi_1_3_glnxa64');
%addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/poly_helpers');
addpath('/home/research/g.dehghanpoor/libraries/libraries_in_usage/SOSTOOLS-3.300/internal');
addpath('/home/research/g.dehghanpoor/libraries/sdpa-c/mex')
addpath('/home/research/g.dehghanpoor/libraries/sdpa-c');
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mfiles')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mfiles/writeFunctions')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mex')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V200')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V210SubPrograms')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/V260SubPrograms')
addpath('/home/research/g.dehghanpoor/libraries/SparsePOP300/subPrograms/Mfiles/Solvers')
folder = '/home/research/g.dehghanpoor/libraries/Sparse_BSOS-master';
genpath(folder)
addpath(genpath(folder));


begin=0;


programRelaxOrder=1;
for i=1:20
    multiple_weights_fixing;


     pos=0;
    if hd<5
        pos=1;
    end

POP=0;
    add_negative(begin,pos,POP,mon_num);
    begin=0;


%     
end



