%% load raw data
clear filename;
datapath = ['E:\Brandeis\data\20200612\results\'];
testpath = ['E:\Brandeis\data\annotation_b202006\'];
id = 1; gap = []; minus2H = [];
filename{id} = 'O_4.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_5.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_6.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_7.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_8.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_9.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_10.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_11.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_12.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_13.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_14.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'O_15.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1


%%
x = 28;
specSet = CSpectrum.empty(length(12), 0);
testSet = CSpectrum.empty(length(12), 0);
idx = 1;
for f = 1 : 12
    load([datapath, replace(filename{f}, '.txt', '.mat')]);
    specSet{idx} = aSpec;
    testSet{idx} = CSpectrum.load_ann([testpath, '20200612_', replace(filename{f}, '.txt', '.ann')]);
    idx = idx + 1;
end
%%
score = lc.rank_candidates(specSet);
cMap = lc.mClassifierMap;
%%
[accurateNum, accuracy, wrongList, wrongData] = lc.testScoreMap(score, specSet, testSet);

