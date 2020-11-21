%% 
massBound = [CMass.H*1.5, 30]; 
massAccuracy = 0.005;
masses = [];
use_simple_method = 0;
mBoostNum = 15;
mHoldOut = 0.2;
savePath = 'E:\Brandeis\data\annotation\2020.11.7_large.txt';
saveMapPath = 'E:\Brandeis\data\annotation\2020.11.7_large_map.txt';
saveTrainResultPath = 'E:\Brandeis\data\annotation\train_result.txt';
%% This dataset use SNAP peak picking

clear filename;
dataset = {'annotation_b202006'};
files = [];
numFiles = 0;
for f = 1 : length(dataset)
    if ismac
        datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', dataset{f}, filesep];
    else
        datapath = ['E:\Brandeis\data\', dataset{f}, filesep];
    end
        files = dir( [datapath, '*.ann'] );
        numFiles = numFiles + length(files);
end
%%
specSet = CSpectrum.empty(length(numFiles), 0);
total = 0;
pass = 0;
for f = 1 : length(dataset)
    if ismac
        datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', dataset{f}, filesep];
    else
        datapath = ['E:\Brandeis\data\', dataset{f}, filesep];
    end
        files = dir( [datapath, '*.ann'] );
        for ff = 1 : length( files )
            filename = [datapath, files(ff).name];
            disp( ['Loading ', filename] );
            specSet{ff} = CSpectrum.load_ann( filename );
            total = total + length(specSet{ff}.mPeaks);
        end
end
disp('Done loading raw data ...');

%% calculate
map = containers.Map('KeyType','char','ValueType','int8');
for s = 1 : length(specSet)
    spec = specSet{s};
    for p = 1 : length(spec.mPeaks)
        peak = spec.mPeaks(p);
        s = [peak.type, ' ',peak.RE, '-', peak.NRE];
        if ~isKey(map, s)
            map(s) = 1;
        else
            map(s) = map(s) + 1;
        end
    end
end
%% merge mass
lc = linkageClassifier;
hp=waitbar(0,'merge masses progress');
for s = 1 : length(specSet)
    spec = specSet{s};
    masses = linkageClassifier.calculate_masses(spec, masses, massBound);
    waitbar(s/length(specSet), hp, ['Calculating context spec ' num2str(s) ' ' num2str(100*s/length(specSet)) '%'])
end
delete(hp);
massFeature = lc.calculate_mass_feature(specSet, masses, massBound, massAccuracy);
%% output
%save(specSet)?
fid = fopen(savePath, 'w');
traindata = [];
fprintf(fid, 'm/z\tz\tintensity\tpeaktype\tion-type\tRE\tNRE\tlinkage\tion-formula\tfeature\tfeatureintensity\n');
for s = 1 : length(specSet)
    spec = specSet{s};
    for p = 1:length(spec.mPeaks)
        peak = spec.mPeaks(p);
        com = '';
        if peak.mIsComplement
            com = 'com';
        end
        if strcmp(peak.type, 'B') == 1 && (strcmp(peak.RE, 'Hex') == 1 && strcmp(peak.NRE, 'HexNAc') == 1 || strcmp(peak.RE, 'Hex') == 1 && strcmp(peak.NRE, 'Hex') == 1 || strcmp(peak.RE, 'HexNAc') == 1 && strcmp(peak.NRE, 'Hex') == 1)
            traindata = [traindata, peak];
            fprintf(fid, '%-f\t%-d\t%-f\t%-s\t%-s\t%-6s\t%-6s\t%-d\t%-15s\t%-s\t%-s\n', peak.mRawMZ, peak.mRawZ, peak.mZscore, com, peak.type, peak.RE, peak.NRE, peak.linkage, peak.mComment, vec2str(peak.mFeature), vec2str(peak.mFeatureIntensity));
        end
    end
end
fclose(fid);
%% output map
fid2 = fopen(saveMapPath, 'w');
keyss = keys(map);
for k = 1 : length(keyss)
    key = keyss{k};
    fprintf(fid2, '%s\t%d\n', key, map(key));
end
fclose(fid2);
%% train linkage classifier
[dataVectors, dataSources] = linkageClassifier.prepare_training_data(traindata);

%% test linkage classifier
vs = nchoosek(keys(dataVectors),2);
fid3 = fopen(saveTrainResultPath, 'w');
fprintf(fid3, '%13s\t%13s\t%13s\t%13s\t%13s\n',"CLASSIFIER","POS","NEG","TRAINERR","TESTERR");
classifierMap = containers.Map;
for m = trainMethods.methodList
    for v = 1 : length(vs)
    pos = vs{v,1};
    neg = vs{v,2};
    posX = dataVectors(pos);
    cutPos = round(size(posX, 1)*0.8);
    posTestX = posX(cutPos:end,:);
    posX = posX(1:(cutPos-1), :);
    posSource = dataSources(pos);
    negX = dataVectors(neg);
    cutNeg = round(size(negX, 1)*0.8);
    negTestX = negX(cutNeg:end,:);
    negX = negX(1:(cutNeg-1), :);
    negSource = dataSources(neg);
    X = [posX; negX];
    Y = [ones(size(posX,1), 1); zeros(size(negX,1), 1)];
    testX = [posTestX; negTestX];
    testY = [ones(size(posTestX,1),1); zeros(size(negTestX,1),1)];
    XInfo = [posSource, negSource];
    weights = [ones(size(posX,1),1)/size(posX,1); ones(size(negX,1),1)/size(negX,1)] / 2;
    hp = cell(2,1);
    hp{1} = mBoostNum;
    hp{2} = weights;

        currclassifier = trainMethods.train_by_method(X, Y, m, hp);
        classifierMap(strcat([pos, neg],m)) = currclassifier;
        trainLabel = predict(currclassifier.classifier, X);
        testLabel = predict(currclassifier.classifier, testX);
        if strcmp(m, 'randomForest') == 1
            trainLabel = str2double(trainLabel);
            testLabel = str2double(testLabel);
        end
        trainError = sum(abs(trainLabel - Y)) / length(Y);
        testError = sum(abs(testLabel - testY)) / length(testY);
        fprintf(fid3, '%13s\t%13s\t%13s\t%13f\t%13f\n', m, pos, neg, trainError, testError);
    end
end
fclose(fid3);
disp('Done!')
%%
save lc;