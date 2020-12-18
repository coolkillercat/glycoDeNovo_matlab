%% load raw data
clear filename;
datapath = ['E:\Brandeis\data\20170615.2\results\'];
testpath = ['E:\Brandeis\data\annotation_b202006\'];
id = 1; gap = []; minus2H = [];
filename{id} = 'Lewis B.O18.Na.EED.PM.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 1
filename{id} = 'Lewis B.O18.Cs.EED.PM.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 2
filename{id} = 'Lewis Y.O18.Na.EED.PM.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 3
filename{id} = 'Lewis Y.O18.Cs.EED.PM.txt'; gap(id) = 1; minus2H(id) = 0; id = id + 1; % 4
filename{id} = 'LNFP 1.O18.Na.EED.PM.txt'; gap(id) = 1; minus2H(id) = 1; id = id + 1; % 5
filename{id} = 'LNFP 1.O18.Cs.EED.PM.txt'; gap(id) = 1; minus2H(id) = 1; id = id + 1; % 6
filename{id} = 'LNFP 2.O18.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 7
filename{id} = 'LNFP 2.O18.Cs.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 8
filename{id} = 'LNFP 3.O18.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 9
filename{id} = 'LNFP 3.O18.Cs.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 10
filename{id} = 'Sialyl Lewis A.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 11
filename{id} = 'Sialyl Lewis X.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 12
filename{id} = 'CelHex.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 13
filename{id} = 'MalHex.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 14
filename{id} = 'LNnT.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 15
filename{id} = 'LNT.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 16
filename{id} = 'Lewis B.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 17
filename{id} = 'Man9.O18.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 18
filename{id} = 'N002.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 19
filename{id} = 'N003.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 20
filename{id} = 'N012.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 21
filename{id} = 'N013.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 22
filename{id} = 'N222.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 23
filename{id} = 'N223.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 24
filename{id} = 'N233.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 25
filename{id} = 'NA2F.O18.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 26
filename{id} = 'LNFP 2.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 27
filename{id} = 'A2F.DR.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 1; id = id + 1; % 28
filename{id} = 'A2F.Reduced.Na.EED.PM.txt'; gap(id) = 0; minus2H(id) = 0; id = id + 1; % 29

%%
x = 28;
specSet = CSpectrum.empty(length(1), 0);
testSet = CSpectrum.empty(length(1), 0);
idx = 1;
for f = 28 : 28
    load([datapath, 'rec.',replace(filename{f}, '.txt', '.mat')]);
    specSet{idx} = specU;
    testSet{idx} = CSpectrum.load_ann([testpath, '20170615.2_', replace(filename{f}, '.txt', '.ann')]);
    idx = idx + 1;
end
%%
score = lc.rank_candidates(specSet);
cMap = lc.mClassifierMap;
%%
[accurateNum, accuracy, wrongList, wrongData] = lc.testScoreMap(score, specSet, testSet);
