%% 2019/08/20
%% skip training flag
flag = 1;
%% IonClassifier -- load training data
if flag == 0
spectraTrain = load_saved_spectra( '20170615.2' ); % or use '20170615'

%% IonClassifier -- train classifier
ionClassifier = CIonClassifier;
ionClassifier.mMassAccuracy = 0.01;
[trainVectors, trainData] = ionClassifier.train( spectraTrain );
disp( 'Train IonClassifier done.' );
end
%% save ion classifier here

fid = fopen( 'ionclassifier_massfeatures.txt', 'w');
fprintf(fid, '%.15f\n', ionClassifier.mMassFeatures );
fclose(fid);


%% save classifier here

posSet = {'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C'};
negSet = {'C', 'Y', 'Z', 'O', 'B', 'Y', 'Z', 'O'};
for round = 1 : 8
    posIon = posSet{round};
    negIon = negSet{round};
    p_v_n = [posIon, '_v_', negIon];
    fid = fopen( ['ionclassifier_', p_v_n, '.txt'], 'w');
    for k = 1 : 100
        fprintf(fid, 'Model %d\n', k);
        fprintf(fid, '\tNumNodes %d\n', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.Trained{k}.NumNodes);
        fprintf(fid, '\tChildren\n');
        fprintf(fid, '\t\t%d %d\n', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.Trained{k}.Impl.Children');
        fprintf(fid, '\tEndChildren\n');
        fprintf(fid, '\tClassProb\n');
        fprintf(fid, '\t\t%.15f %.15f\n', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.Trained{k}.Impl.ClassProb');
        fprintf(fid, '\tEndClassProb\n');
        fprintf(fid, '\tCutPoint\n\t\t');
        fprintf(fid, '%.15f ', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.Trained{k}.Impl.CutPoint);
        fprintf(fid, '\n\tCutVar\n\t\t');
        fprintf(fid, '%d ', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.Trained{k}.Impl.CutVar);
        fprintf(fid, '\n\n');
    end
    fprintf(fid, 'TrainedWeights\n');
    fprintf(fid, '%.15f ', ionClassifier.mClassifier.(p_v_n).classifier.Trained{1}.TrainedWeights);
    fclose(fid);
end
%% Remove unused features
featureIdxes = [];
clfs = fields( ionClassifier.mClassifier );
for k = 1 : length(clfs)
    featureIdxes = [featureIdxes, find(ionClassifier.mClassifier.(clfs{k}).FeatureImportance)];
end
numMassFeatures = length(ionClassifier.mMassFeatures);
idxes = find( featureIdxes > numMassFeatures );
featureIdxes(idxes) = featureIdxes(idxes) - numMassFeatures;
featureIdxes = unique(featureIdxes);
idxes = find( featureIdxes <= numMassFeatures );
featureIdxes = featureIdxes(idxes);
massFeatures = ionClassifier.mMassFeatures(featureIdxes);

%% Check 
% ionClassifier.mClassifier.B_v_Y.classifier.Trained{1}.TrainedWeights
% ionClassifier.mClassifier.B_v_Y.classifier.Trained{1}.Trained{k}
% edit classreg.learning.classif.ClassificationModel
% ionClassifier.mClassifier.B_v_Y.classifier.Trained{1}.Impl

%% classificationtreescore
testData = ionClassifier.mTestData;
testVectors = ionClassifier.mTestVectors;
fid = fopen('treescore', 'w');
fid2 = fopen('save_testVectors.txt', 'w');
for ion = 'BC'
                if ion == 'B'
                    vs_BC = [ion, '_v_C'];
                else
                    vs_BC = [ion, '_v_B'];
                end
                vs_Y = [ion, '_v_Y'];
                vs_Z = [ion, '_v_Z'];
                vs_O = [ion, '_v_O'];
                
                for k = 1 : length(testData.(ion))
                    filename = ['matlab_', ion, k, '.txt'];
                    fid3 = fopen(['matlab_', ion, num2str(k), '.txt'], 'w');
                    fprintf(fid2, '%c%d\n', ion, k);
                    fprintf(fid2, '%.16f\n', testVectors.(ion)(k,:));
                    fprintf(fid2, '\n\n\n');
                    fprintf(fid3, '%.16f\n', testVectors.(ion)(k,:));
                    fclose(fid3);
                    key = uint32(testData.(ion)(k).spectrumID * 10000 + testData.(ion)(k).peakID);
                    for m = 1 : 100
                    % Calculate its score of being a B-ion
                    ws = zeros(1, 4);
                    [~, b] = ionClassifier.mClassifier.(vs_BC).classifier.Trained{1}.Trained{m}.predict( testVectors.(ion)(k,:) );
                    ws(1) = b(2);
                    if strcmp( ionClassifier.mClassifier.(vs_BC).classifier.CrossValidatedModel, 'Bag' )
                        ws(1) = ws(1) * ionClassifier.mClassifier.(vs_BC).classifier.NumTrainedPerFold;
                    end
                    
                    [~, b] = ionClassifier.mClassifier.(vs_Y).classifier.Trained{1}.Trained{m}.predict( testVectors.(ion)(k,:) );
                    ws(2) = b(2);
                    if strcmp( ionClassifier.mClassifier.(vs_Y).classifier.CrossValidatedModel, 'Bag' )
                        ws(2) = ws(2) * ionClassifier.mClassifier.(vs_Y).classifier.NumTrainedPerFold;
                    end
                    [~, b] = ionClassifier.mClassifier.(vs_Z).classifier.Trained{1}.Trained{m}.predict( testVectors.(ion)(k,:) );
                    ws(3) = b(2);
                    if strcmp( ionClassifier.mClassifier.(vs_Z).classifier.CrossValidatedModel, 'Bag' )
                        ws(3) = ws(3) * ionClassifier.mClassifier.(vs_Z).classifier.NumTrainedPerFold;
                    end
                    [~, b] = ionClassifier.mClassifier.(vs_O).classifier.Trained{1}.Trained{m}.predict( testVectors.(ion)(k,:) );
                    ws(4) = b(2);
                    if strcmp( ionClassifier.mClassifier.(vs_O).classifier.CrossValidatedModel, 'Bag' )
                        ws(4) = ws(4) * ionClassifier.mClassifier.(vs_O).classifier.NumTrainedPerFold;
                    end
                    %ionScoreMap.(ion)(key) = min(ws);
                    fprintf(fid, '\nws: %f %f %f %f\nion: %c model: %d key: %d k: %d\nscore: %f\n\n', ws(1), ws(2), ws(3), ws(4), ion, m, key, k, min(ws));
                    fprintf(fid, '%s\n', vs_BC);
                    fprintf(fid, 'cutpoint: \t%f %f %f\n', ionClassifier.mClassifier.(vs_BC).classifier.Trained{1}.Trained{m}.Impl.CutPoint);
                    fprintf(fid, 'cutvar: \t%d %d %d\n', ionClassifier.mClassifier.(vs_BC).classifier.Trained{1}.Trained{m}.Impl.CutVar);
                    fprintf(fid, 'value: \t%f\n', testVectors.(ion)(k, ionClassifier.mClassifier.(vs_BC).classifier.Trained{1}.Trained{m}.Impl.CutVar(1)));
                    fprintf(fid, 'prob: \t%f %f %f\n\t%f %f %f', ionClassifier.mClassifier.(vs_BC).classifier.Trained{1}.Trained{m}.Impl.ClassProb);
                    fprintf(fid, '\n%s\n', vs_Y);
                    fprintf(fid, 'cutpoint: \t%f %f %f\n', ionClassifier.mClassifier.(vs_Y).classifier.Trained{1}.Trained{m}.Impl.CutPoint);
                    fprintf(fid, 'cutvar: \t%d %d %d\n', ionClassifier.mClassifier.(vs_Y).classifier.Trained{1}.Trained{m}.Impl.CutVar);
                    fprintf(fid, 'value: \t%f\n', testVectors.(ion)(k, ionClassifier.mClassifier.(vs_Y).classifier.Trained{1}.Trained{m}.Impl.CutVar(1)));
                    fprintf(fid, 'prob: \t%f %f %f\n\t%f %f %f', ionClassifier.mClassifier.(vs_Y).classifier.Trained{1}.Trained{m}.Impl.ClassProb);
                    fprintf(fid, '\n%s\n', vs_Z);
                    fprintf(fid, 'cutpoint: \t%f %f %f\n', ionClassifier.mClassifier.(vs_Z).classifier.Trained{1}.Trained{m}.Impl.CutPoint);
                    fprintf(fid, 'cutvar: \t%d %d %d\n', ionClassifier.mClassifier.(vs_Z).classifier.Trained{1}.Trained{m}.Impl.CutVar);
                    fprintf(fid, 'value: \t%f\n', testVectors.(ion)(k, ionClassifier.mClassifier.(vs_Z).classifier.Trained{1}.Trained{m}.Impl.CutVar(1)));
                    fprintf(fid, 'prob: \t%f %f %f\n\t%f %f %f', ionClassifier.mClassifier.(vs_Z).classifier.Trained{1}.Trained{m}.Impl.ClassProb);
                    fprintf(fid, '\n%s\n', vs_O);
                    fprintf(fid, 'cutpoint: \t%f %f %f\n', ionClassifier.mClassifier.(vs_O).classifier.Trained{1}.Trained{m}.Impl.CutPoint);
                    fprintf(fid, 'cutvar: \t%d %d %d\n', ionClassifier.mClassifier.(vs_O).classifier.Trained{1}.Trained{m}.Impl.CutVar);
                    fprintf(fid, 'value: \t%f\n', testVectors.(ion)(k, ionClassifier.mClassifier.(vs_O).classifier.Trained{1}.Trained{m}.Impl.CutVar(1)));
                    fprintf(fid, 'prob: \t%f %f %f\n\t%f %f %f', ionClassifier.mClassifier.(vs_O).classifier.Trained{1}.Trained{m}.Impl.ClassProb);
                    end
                  
                end
end
fclose(fid);
fclose(fid2);

