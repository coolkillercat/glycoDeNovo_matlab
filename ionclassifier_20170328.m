%%
% Use data with reducing end moficiations to train IonClassifier, and apply
% the trained IonClassifier to the data without reducing end moficiation.
%%
spectra = load_saved_spectra( );

%%
flagTrain = zeros(1, length(spectra));
flagTest = zeros(1, length(spectra));
for k = 1 : length(flagTrain)
    flagTrain(k) = strcmp( spectra(k).mExperimentMethod, 'EED' ) & ~strcmp( spectra(k).mReducingEndModification, '' ) & ~strcmp( spectra(k).mMetal, 'Proton' );
    flagTest(k) = strcmp( spectra(k).mExperimentMethod, 'EED' ) & strcmp( spectra(k).mReducingEndModification, '' ) & ~strcmp( spectra(k).mMetal, 'Proton' );
end
spectraTrain = spectra(flagTrain > 0);
spectraTest = spectra(flagTest > 0);

%%
[vectors, data] = prepare_training_data( spectra, 0.02, 0, 0 );

%%
trainData = data; trainVectors = vectors;
testData = data; testVectors = vectors;
ions = {'B', 'C', 'Y', 'Z', 'O'};
for k = 1 : length(ions)
    ion = ions{k};
    disp(ion);
    flagA = ones(1, length(data.(ion)));
    flagB = ones(1, length(data.(ion)));
    for m = 1 : length(data.(ion))
        flagA(m) = flagTrain( data.(ion)(m).glycanID );
        flagB(m) = flagTest( data.(ion)(m).glycanID );
    end
    trainData.(ion) = data.(ion)(flagA > 0);
    testData.(ion) = data.(ion)(flagB > 0);
    trainVectors.(ion) = vectors.(ion)(flagA > 0, :);
    testVectors.(ion) = vectors.(ion)(flagB > 0, :);
end

%%
ionClassifier = train_ionclassifier( trainVectors, trainData, 0 );
results = predict_by_ionclassifier( spectraTest, testData, testVectors, ionClassifier );

%% Don't need the codes below

%%
posSet = {'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C', 'Y', 'Y', 'Y', 'Y', 'Z', 'Z', 'Z', 'Z'}; % 16
negSet = {'C', 'Y', 'Z', 'O', 'B', 'Y', 'Z', 'O', 'B', 'C', 'Z', 'O', 'C', 'B', 'Y', 'O'};
do_LOO = 0;
LOOResult = [];
for round = 1 : length(posSet)

    posIon = posSet{round};
    negIon = negSet{round}; % no enough data, put them in one bucket.
    % use_simple_model = {'B_v_Y', 'B_v_Z', 'C_v_Y'};
    use_simple_model = {};

    p_v_n = [posIon, '_v_', negIon];
    disp( p_v_n );
    
    posX = []; negX = [];
    
    % find samples that appear in both posX and negX
    posSource = [];
    negSource = [];
    
    for n = posIon
        posX = [posX; trainVectors.(n)];
        posSource = [posSource; [[trainData.(n).glycanID]; [trainData.(n).peakID]]'];
    end
    
    for n = negIon
        negX = [negX; trainVectors.(n)];
        negSource = [negSource; [[trainData.(n).glycanID]; [trainData.(n).peakID]]'];
    end
    
    X = [posX; negX];
    Y = [ones(size(posX,1), 1); zeros(size(negX,1), 1)];
    XInfo = [posSource; negSource];
    weights = [ones(size(posX,1),1)/size(posX,1); ones(size(negX,1),1)/size(negX,1)] / 2;
    
    if do_LOO
        % Leave-one-out
        N = size(X,1);
        LOO = NaN(2, N);
        accP = 0; accN = 0; numP = 0; numN = 0;
        boostNum = 100;
        holdOut = 0.20;
        
        for k = 1 : N
            flag = ones(1, N);
            % flag(k) = 0;
            flag((XInfo(k,1) == XInfo(:,1)) & (XInfo(k,2) == XInfo(:,2))) = 0; % some peaks can be mutliple types (B, C, Y, Z)
            flag = flag > 0;
            if any(strcmp( p_v_n, use_simple_model ) )
                tree = fitctree(X(flag,:),Y(flag));
                [a, b] = tree.predict(X(k,:));
                LOO(1, k) = a;
                LOO(2, k) = b(a+1);
            else
                ClassTreeEns = fitensemble(X(flag,:), Y(flag), 'AdaBoostM1', boostNum, 'Tree', 'Holdout', holdOut, 'Weights', weights(flag));
                if isempty(ClassTreeEns.Trained{1}.predictorImportance)
                    tree = fitctree(X, Y);
                    [a, b] = tree.predict(X(k,:));
                    LOO(1, k) = a;
                    LOO(2, k) = b(a+1);
                else
                    [a, b] = ClassTreeEns.Trained{1}.predict(X(k,:));
                    LOO(1, k) = a;
                    LOO(2, k) = b(a+1);
                end
            end
            if Y(k) == 1
                accP = accP + (Y(k) == LOO(1, k));
                numP = numP + 1;
            else
                accN = accN + (Y(k) == LOO(1, k));
                numN = numN + 1;
            end
            msg = sprintf( '%s vs %s: (%4d,%4d)\t[%d --> %d (%5.2f)]\t%d/%d\t%d/%d\n', posIon, negIon, ...
                XInfo(k,1), XInfo(k,2), Y(k), LOO(1,k), LOO(2,k), accP, numP, accN, numN );
            fprintf( msg );
        end
        disp( [p_v_n, ' LOO finished.'] );
    end
    
    temp = [trainVectors.massFeatures, trainVectors.massFeatures, trainVectors.massFeatures];
    LOOResult.(p_v_n).X = X;
    LOOResult.(p_v_n).Y = Y;
    LOOResult.(p_v_n).XInfo = XInfo;
    LOOResult.(p_v_n).BoostingNum = boostNum;
    LOOResult.(p_v_n).BoostingHoldout = holdOut;
    if do_LOO
        LOOResult.(p_v_n).LOO = LOO;
        LOOResult.(p_v_n).fail = XInfo( LOO(1,:)' ~= Y, : );
    end
    
    % remove shared samples
    [~, ia] = unique( XInfo, 'rows' );
    
    if any(strcmp( p_v_n, use_simple_model ) )
        LOOResult.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
        LOOResult.(p_v_n).FeatureImportance = LOOResult.(p_v_n).classifier.predictorImportance;
    else
        LOOResult.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'AdaBoostM1', boostNum, 'Tree', 'Holdout', holdOut, 'Weights', weights(ia));
        if isempty( LOOResult.(p_v_n).classifier.Trained{1}.predictorImportance )
            LOOResult.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
            LOOResult.(p_v_n).FeatureImportance = LOOResult.(p_v_n).classifier.predictorImportance;
        else
            LOOResult.(p_v_n).FeatureImportance = LOOResult.(p_v_n).classifier.Trained{1}.predictorImportance;
        end
    end
    LOOResult.(p_v_n).SelectedFeatureIdx = find( LOOResult.(p_v_n).FeatureImportance(1:end-2) > 0 );
    LOOResult.(p_v_n).SelectedFeatureMass = temp( LOOResult.(p_v_n).SelectedFeatureIdx );
    disp( LOOResult.(p_v_n) );
    % disp( num2cell([LOOResult.(p_v_n).SelectedFeatureIdx; LOOResult.(p_v_n).SelectedFeatureMass]) );
    
end % round

% save LOOResult_20170327 LOOResult trainVectors trainData;

%%
num = length( spectra );
ICScores = cell(1, num);
CVResults = {'SpecID', 'Name', 'Formula', 'Fragmentation', 'Precursor', 'ReducingEndModification', 'Metal', '# Peaks', '# Interpretable Peaks', '# Reconstructed Peaks', '# Candidates', 'Rank by Supporting Peak Number', 'Rank by IonClassifier', 'File'};
for s = 1 : num
    if isempty(spectra(s).mPeaks(end).mInferredSuperSet)
        continue;
    end
    if strcmp(spectra(s).mExperimentMethod, 'EED' ) == 0
        continue;
    end

    TSS = spectra(s).mPeaks(end).mInferredSuperSet;
    ICScores{s} = cell(1, length(TSS.mTopologies));
    for t = 1 : length(TSS.mTopologies)
        tp = TSS.mTopologies(t);
        ICScores{s}{t} = zeros(4, length(tp.mSupportPeaks)-1);
        weights = zeros(1, length(tp.mSupportPeaks)-1);
        for m = 1 : length(tp.mSupportPeaks(1:end-1))
            peakID = tp.mSupportPeaks(m);
            ICScores{s}{t}(1, m) = peakID;
            type = 0;
            for peakTSS = spectra(s).mPeaks(peakID).mInferredSuperSet
                if isempty( peakTSS.mFormulas ), continue; end
                switch peakTSS.mTargetPeaks(2, peakTSS.mTargetPeaks(1, :) == peakID)
                    case {1, 11, 21}
                        if spectra(s).mPeaks(peakID).mComplement < 0
                            ICScores{s}{t}(2, m) = -spectra(s).mPeaks(peakID).mComplement;
                            type = -1;
                        else
                            type = 1;
                        end
                    case {2, 12, 22}
                        if spectra(s).mPeaks(peakID).mComplement < 0
                            ICScores{s}{t}(2, m) = -spectra(s).mPeaks(peakID).mComplement;
                            type = -2;
                        else
                            type = 2;
                        end
                end
            end
            
            if type == 1
                w = zeros(1, 4);
                missed = 0;
                idx = find( (LOOResult.B_v_Y.XInfo(:,1) == s) & (LOOResult.B_v_Y.XInfo(:,2) == peakID) );
                if isempty( idx )
                    missed = missed + 1;
                    w(1) = 0;
                else
                    % w(1) = LOOResult.B_v_Y.LOO(2, idx);
                    w(1) = (2 * LOOResult.B_v_Y.LOO(1, idx) - 1) * LOOResult.B_v_Y.LOO(2, idx);
                end
                idx = find( (LOOResult.B_v_Z.XInfo(:,1) == s) & (LOOResult.B_v_Z.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.B_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.B_v_Z.LOO(1, idx) - 1) * LOOResult.B_v_Z.LOO(2, idx);
                end
                idx = find( (LOOResult.B_v_C.XInfo(:,1) == s) & (LOOResult.B_v_C.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.B_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.B_v_C.LOO(1, idx) - 1) * LOOResult.B_v_C.LOO(2, idx);
                end
                idx = find( (LOOResult.B_v_O.XInfo(:,1) == s) & (LOOResult.B_v_O.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.B_v_O.LOO(2, idx);
                    w(4) = (2 * LOOResult.B_v_O.LOO(1, idx) - 1) * LOOResult.B_v_O.LOO(2, idx);
                end
                if missed == 4
                    fprintf( 'B Missed [%d %d]\n', s, peakID );
                end
                weights(m) = min(w);
            elseif type == -1
                w = [0, 0, 0, 0];
                missed = 0;
                idx = find( (LOOResult.Y_v_O.XInfo(:,1) == s) & (LOOResult.Y_v_O.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty( idx )
                    missed = missed + 1;
                    w(1) = 0;
                else
                    % w(1) = LOOResult.Y_v_O.LOO(2, idx);
                    w(1) = (2 * LOOResult.Y_v_O.LOO(1, idx) - 1) * LOOResult.Y_v_O.LOO(2, idx);
                end
                idx = find( (LOOResult.Y_v_Z.XInfo(:,1) == s) & (LOOResult.Y_v_Z.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.Y_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.Y_v_Z.LOO(1, idx) - 1) * LOOResult.Y_v_Z.LOO(2, idx);
                end
                idx = find( (LOOResult.Y_v_C.XInfo(:,1) == s) & (LOOResult.Y_v_C.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.Y_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.Y_v_C.LOO(1, idx) - 1) * LOOResult.Y_v_C.LOO(2, idx);
                end
                idx = find( (LOOResult.Y_v_B.XInfo(:,1) == s) & (LOOResult.Y_v_B.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.Y_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.Y_v_B.LOO(1, idx) - 1) * LOOResult.Y_v_B.LOO(2, idx);
                end
                weights(m) = min(w);
                if missed == 4
                    fprintf( 'Y Missed [%d %d]\n', s, peakID );
                end
            elseif type == 2
                w = [0 0 0 0];
                missed = 0;
                idx = find( (LOOResult.C_v_Y.XInfo(:,1) == s) & (LOOResult.C_v_Y.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(1) = 0;
                else
                    % w(1) = LOOResult.C_v_Y.LOO(2, idx);
                    w(1) = (2 * LOOResult.C_v_Y.LOO(1, idx) - 1) * LOOResult.C_v_Y.LOO(2, idx);
                end
                idx = find( (LOOResult.C_v_Z.XInfo(:,1) == s) & (LOOResult.C_v_Z.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.C_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.C_v_Z.LOO(1, idx) - 1) * LOOResult.C_v_Z.LOO(2, idx);
                end
                idx = find( (LOOResult.C_v_O.XInfo(:,1) == s) & (LOOResult.C_v_O.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.C_v_O.LOO(2, idx);
                    w(3) = (2 * LOOResult.C_v_O.LOO(1, idx) - 1) * LOOResult.C_v_O.LOO(2, idx);
                end
                idx = find( (LOOResult.C_v_B.XInfo(:,1) == s) & (LOOResult.C_v_B.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.C_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.C_v_B.LOO(1, idx) - 1) * LOOResult.C_v_B.LOO(2, idx);
                end
                
                if missed == 4
                    fprintf( 'C Missed [%d %d]\n', s, peakID );
                end
                weights(m) = min(w);
            elseif type == -2
                w = [0 0 0 0];
                missed = 0;
                idx = find( (LOOResult.Z_v_O.XInfo(:,1) == s) & (LOOResult.Z_v_O.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty( idx )
                    missed = missed + 1;
                    w(1) = 0;
                else
                    % w(1) = LOOResult.Z_v_O.LOO(2, idx);
                    w(1) = (2 * LOOResult.Z_v_O.LOO(1, idx) - 1) * LOOResult.Z_v_O.LOO(2, idx);
                end
                idx = find( (LOOResult.Z_v_Y.XInfo(:,1) == s) & (LOOResult.Z_v_Y.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.Z_v_Y.LOO(2, idx);
                    w(2) = (2 * LOOResult.Z_v_Y.LOO(1, idx) - 1) * LOOResult.Z_v_Y.LOO(2, idx);
                end
                idx = find( (LOOResult.Z_v_C.XInfo(:,1) == s) & (LOOResult.Z_v_C.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.Z_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.Z_v_C.LOO(1, idx) - 1) * LOOResult.Z_v_C.LOO(2, idx);
                end
                idx = find( (LOOResult.Z_v_B.XInfo(:,1) == s) & (LOOResult.Z_v_B.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.Z_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.Z_v_B.LOO(1, idx) - 1) * LOOResult.Z_v_B.LOO(2, idx);
                end
                weights(m) = min(w);
                if missed == 4
                    fprintf( 'Z Missed [%d %d]\n', s, peakID );
                end
            end
             ICScores{s}{t}(3, m) = weights(m);
             ICScores{s}{t}(4, m) = type;
        end
        tp.mScore = sum(weights);
        % tp.mScore = sum( 2 * (weights-0.5) ) + 1;
    end
    
    numInterpretable = 0;
    numReconstructed = 0;
    for m = 1 : length(spectra(s).mPeaks)
        if ~isempty(spectra(s).mPeaks(m).mInferredSuperSet)
            numInterpretable = numInterpretable + 1;
            if ~isempty(spectra(s).mPeaks(m).mInferredFormulas)
                numReconstructed = numReconstructed + 1;
            end
        end
    end
    if numReconstructed == 0, continue; end
    idx = strfind( spectra(s).comment, '.' );
    if ~isempty(idx)
        formula = spectra(s).comment(idx+2:end);
    else
        formula = spectra(s).comment;
    end
    formula = strtrim( formula );
    formula = strrep( formula, '][', '] [');
    
    name = spectra(s).filename;
    idx = strfind( name, '\' );
    name = spectra(s).filename(idx(end)+1:end);
    idx = strfind( name, ' ' );
    if ~isempty( idx )
        name = name(5:idx(1)-1);
    else
        name = strrep( name(5:end), '.mat', '' );
    end
    
    % Get groundtruth
    g = CGlycan;  g.parse( formula );   g.order_branch_by_name();
    
    % get rank information
    rankInfo = spectra(s).get_rank( g.inferred_formula );
    
    if rankInfo.rankBySupportingPeak == 0, continue; end
    
    numPeaks = [num2str(length(spectra(s).mPeaks)), '(', num2str(sum( [spectra(s).mPeaks.mComplement] == 0)), ')'];
    byPeaks = [num2str(rankInfo.rankBySupportingPeak), '(', num2str(rankInfo.corankBySupportingPeak-1), ')'];
    byIonClassifier = [num2str(rankInfo.rankByIonclassifier), '(', num2str(rankInfo.corankByIonclassifier-1), ')'];
    disp( [byPeaks, ', ', byIonClassifier, ', ', formula, ': ', spectra(s).filename] );
    CVResults(end+1,:) = {s, name, formula, spectra(s).mExperimentMethod, spectra(s).mPrecursor, ...
        spectra(s).mReducingEndModification, spectra(s).mMetal, ...
        numPeaks, numInterpretable, numReconstructed, ...
        length(spectra(s).mPeaks(end).mInferredFormulas), ...
        byPeaks, byIonClassifier, spectra(s).filename };
end
disp( 'Finished' );

%%
% [~, idx] = sort( [CVResults{2:end, 5}] );
% CVResults(2:end, :) = CVResults(idx+1,:);
xlswrite( 'LOO_results.EED.xlsx', CVResults );
% save ionclassifier.intensity.EED_O18.mat trainData trainVectors LOOResult ICScores spectra
%% Inspect the selected features
p_v_n = 'C_v_O';
num = length( trainVectors.massFeatures );
sf.idx = LOOResult.(p_v_n).SelectedFeatureIdx;
sf.importance = LOOResult.(p_v_n).FeatureImportance(sf.idx);
sf.pos = LOOResult.(p_v_n).X(LOOResult.(p_v_n).Y > 0, sf.idx);
sf.neg = LOOResult.(p_v_n).X(LOOResult.(p_v_n).Y == 0, sf.idx);
idx = sf.idx; idx(idx > num) = idx(idx>num)-num;
sf.massFeatures = trainVectors.massFeatures(idx);
[sf.importance, idx] = sort( sf.importance, 'descend' );
sf.idx = sf.idx(idx);
sf.massFeatures = sf.massFeatures(idx);
sf.pos = sf.pos(:, idx);
sf.neg = sf.neg(:, idx);

%% Visualize important features
figure;
clear a;
k = 1;
num = length( trainVectors.massFeatures );
for r = 1 : 4
    for c = 1 : 5
        if k > length(sf.massFeatures)
            break;
        end
        a{1} = sf.pos(:,k);
        a{2} = sf.neg(:,k);
        if sf.idx(k) > num && sf.idx(k) < num+num
            a{1} = a{1}( sf.pos(:,k) ~= 0 );
            a{2} = a{2}( sf.neg(:,k) ~= 0 );
        end
        subplot(4, 5, k);
        if sf.idx(k) <= num
            if k == 1
                compare_hist(a, 2, 1, {'P', 'N'});
            else
                compare_hist(a, 2, 1);
            end
        else
            if k == 1
                compare_hist(a, [], 1, {'P', 'N'});
            else
                compare_hist(a, [], 1);
            end
        end
        title( ['(', num2str(sf.massFeatures(k)), ', ', num2str(sf.importance(k)), ')'] );
        k = k + 1;
    end
end

%%

%% LASSO
[coeff, fitInfo] = lasso( X, Y, 'CV', 10, 'Weights', weights);
lassoPlot(coeff,fitInfo,'PlotType','CV');
% lassoPlot(coeff,fitInfo,'PlotType','Lambda','XScale','log');
title( 'LASSO' );

%% Decision tree
test = fitctree(X,Y);
view(test,'Mode','graph')


%% Boosting Tree
ClassTreeEns = fitensemble(X,Y,'AdaBoostM1',1,'Tree', 'Holdout',0.2, 'Weights', weights);
genError = kfoldLoss(ClassTreeEns,'Mode','Cumulative');
figure;
plot(genError);
title('AdaBoost Tree');

