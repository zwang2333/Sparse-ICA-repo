% test the sparsefastica function

%sparse fastICA  sim
%addpath("/Users/ben/Desktop/work2/DICA/2stageNEW/MIToolbox")
%addpath("/Users/zihang/Library/CloudStorage/OneDrive-EmoryUniversity/Research/Ben_RA/SparseICA/Simulations/Sim123/matlab/")

approach = 'symm';
numOfIC = 3;
finetune = 'off';
myy = 1;
stabilization = 'on';
epsilon1 = 0.001;
epsilon2 = 0.000001;
maxNumIterations = 1000;
maxFinetune = 100;
initState = 'rand';
guess = 0;
sampleSize = 1;
displayMode = 'off';
displayInterval = 1;
s_verbose = 'on';
sigma = 1;
rou = 0.965;
func = 'tanh';
round_2nd = 20;

% low SNR situation
load('../../Data/sim123_sparse_lowSNR.mat');

for i = 1:100
    xmat = xmat_matlab_lowSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    [E, D] = pcamat(xmat',1,3);
    [nv, wm, dwm] = whitenv(xmat', E, D);

    X = nv;
    whiteningMatrix = wm;
    dewhiteningMatrix = dwm;

    tStart = tic; 
    [A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, finetune, myy, stabilization, epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, guess, sampleSize, displayMode, displayInterval, s_verbose, sigma, rou, func, round_2nd);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/sparsefast/low/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','A','tEnd','xmat');
end

% medium SNR situation
load('../../Data/sim123_sparse_mediumSNR.mat');

for i = 1:100
    xmat = xmat_matlab_mediumSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    [E, D] = pcamat(xmat',1,3);
    [nv, wm, dwm] = whitenv(xmat', E, D);

    X = nv;
    whiteningMatrix = wm;
    dewhiteningMatrix = dwm;
    
    tStart = tic; 
    [A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, finetune, myy, stabilization, epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, guess, sampleSize, displayMode, displayInterval, s_verbose, sigma, rou, func, round_2nd);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/sparsefast/medium/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','A','tEnd','xmat');
end

% high SNR situation
load('../../Data/sim123_sparse_highSNR.mat');

for i = 1:100
    xmat = xmat_matlab_highSNR((1+1089*(i-1)):(1089+1089*(i-1)),:);
    [E, D] = pcamat(xmat',1,3);
    [nv, wm, dwm] = whitenv(xmat', E, D);

    X = nv;
    whiteningMatrix = wm;
    dewhiteningMatrix = dwm;
    
    tStart = tic; 
    [A, W, time] = sparsefastica(X, whiteningMatrix, dewhiteningMatrix, approach, numOfIC, finetune, myy, stabilization, epsilon1, epsilon2, maxNumIterations, maxFinetune, initState, guess, sampleSize, displayMode, displayInterval, s_verbose, sigma, rou, func, round_2nd);
    tEnd = toc(tStart);
    myS = W * xmat';
    filenm = ['../../Results/sparsefast/high/estS_' num2str(i) '.mat' ];
    save(filenm,'myS','W','A','tEnd','xmat');
end

