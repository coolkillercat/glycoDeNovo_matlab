classdef CIonClassifier < handle  % By Pengyu Hong @ Brandeis University
    properties
        mBoostNum = 100;
        mHoldOut = 0.2;
        
        mMassFeatures = [];
        mMassAccuracy = 0.005;
        
        mClassifier = [];
        mUseOriginalPeaks = 0;
        mUseComplementFlag = 0;
        mEnableX15 = 0;
        
        mTestVectors = [];
        mTestData = [];
    end
    
    methods
        function ionScoreMap = predict_ions(obj, testVectors, testData)
            obj.mTestVectors = testVectors;
            obj.mTestData = testData;
            fid = fopen( ['scoremap.txt'], 'w');
            idxX15_B = find( abs( obj.mMassFeatures - (-27.9949) ) < 0.005 );
            idxX15_C = find( abs( obj.mMassFeatures - (-27.9949 - CMass.H2O) ) < 0.005 );

            ionScoreMap.comment = 'key = uint32(glycanID * 10000 + peakID)';
            ionScoreMap.B = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
            ionScoreMap.C = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
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
                    key = uint32(testData.(ion)(k).spectrumID * 10000 + testData.(ion)(k).peakID);
                    
                    % Calculate its score of being a B-ion
                    ws = zeros(1, 4);
                    [~, b] = obj.mClassifier.(vs_BC).classifier.Trained{1}.predict( testVectors.(ion)(k,:) );
                    ws(1) = b(2);
                    if strcmp( obj.mClassifier.(vs_BC).classifier.CrossValidatedModel, 'Bag' )
                        ws(1) = ws(1) * obj.mClassifier.(vs_BC).classifier.NumTrainedPerFold;
                    end
                    fprintf(fid, '\nb: %f ', b(2));
                    [~, b] = obj.mClassifier.(vs_Y).classifier.Trained{1}.predict( testVectors.(ion)(k,:) );
                    ws(2) = b(2);
                    fprintf(fid, '%f ', b(2));
                    if strcmp( obj.mClassifier.(vs_Y).classifier.CrossValidatedModel, 'Bag' )
                        ws(2) = ws(2) * obj.mClassifier.(vs_Y).classifier.NumTrainedPerFold;
                    end
                    [~, b] = obj.mClassifier.(vs_Z).classifier.Trained{1}.predict( testVectors.(ion)(k,:) );
                    ws(3) = b(2);
                    fprintf(fid, '%f ', b(2));
                    if strcmp( obj.mClassifier.(vs_Z).classifier.CrossValidatedModel, 'Bag' )
                        ws(3) = ws(3) * obj.mClassifier.(vs_Z).classifier.NumTrainedPerFold;
                    end
                    [~, b] = obj.mClassifier.(vs_O).classifier.Trained{1}.predict( testVectors.(ion)(k,:) );
                    ws(4) = b(2);
                    fprintf(fid, '%f \n', b(2));
                    if strcmp( obj.mClassifier.(vs_O).classifier.CrossValidatedModel, 'Bag' )
                        ws(4) = ws(4) * obj.mClassifier.(vs_O).classifier.NumTrainedPerFold;
                    end
                    ionScoreMap.(ion)(key) = min(ws);
                    fprintf(fid, 'ws: %f %f %f %f\nion: %c key: %d\nscore: %f\n\n', ws(1), ws(2), ws(3), ws(4), ion, key, min(ws));
                    % S: 12/30/2018. Test X15
                    if obj.mEnableX15
                        if ion == 'B' && any( testVectors.(ion)(k, idxX15_B) )
                            ionScoreMap.(ion)(key) = ionScoreMap.(ion)(key) + 5;
                        elseif ion == 'C' && any( testVectors.(ion)(k, idxX15_C) )
                            ionScoreMap.(ion)(key) = ionScoreMap.(ion)(key) + 5;
                        end
                    end
                    % E: 12/30/2018. Test X15
                end
            end
            fclose(fid);
        end
        
        function ICScores = rank_candidates(obj, spectra)
            [vectors, signals] = CIonClassifier.extract_data( spectra, obj.mMassFeatures, [], 0.005, obj.mUseOriginalPeaks );
            ionScoreMap = predict_ions(obj, vectors, signals);
            
            for s = 1 : length(spectra)
                if isempty(spectra(s).mPeaks(end).mInferredSuperSet)
                    continue;
                end
                fid = fopen( ['peakweight',num2str(s),'.txt'], 'w');
                TSS = spectra(s).mPeaks(end).mInferredSuperSet;
                ICScores{s} = cell(1, length(TSS.mTopologies));
                for t = 1 : length(TSS.mTopologies)
                    tp = TSS.mTopologies(t);
                    len = length(tp.mSupportPeaks)-1;
                    ICScores{s}{t} = zeros(4, len); % Each peak has [peakID, complementPeakID (if used), type, score]'
                    weights = zeros(1, len);
                    for m = 1 : len
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
                        ICScores{s}{t}(3, m) = type;
                        
                        key = uint32(s * 10000 + peakID);
                        
                        if type == 1 || type == -1
                            weights(m) = ionScoreMap.B(key);
                            spectra(s).mPeaks(peakID).mIonClassifierScores(1) = weights(m);
                        elseif type == 2 || type == -2
                            weights(m) = ionScoreMap.C(key);
                            spectra(s).mPeaks(peakID).mIonClassifierScores(2) = weights(m);
                        end
                        ICScores{s}{t}(4, m) = weights(m);
                        fprintf(fid, '%d\t%d\t%d\t%f\n', ICScores{s}{t}(1,m), ICScores{s}{t}(2,m), ICScores{s}{t}(3,m), ICScores{s}{t}(4,m));
                    end
                    tp.mScore = sum(weights);
                    fprintf(fid, 'tpscore: %f\ntpformula: %s\n\n', tp.mScore, tp.mFormula);
                end
                fclose(fid);
            end
            
            %disp( 'Finished' );
        end
        
        function [trainVectors, trainData] = train(obj, spectra)
            [trainVectors, trainData] = CIonClassifier.prepare_training_data( spectra, [], [], obj.mMassAccuracy, 0, 0 );
            obj.mMassFeatures = trainVectors.massFeatures;
            obj.train_by_ion( trainVectors, trainData );
        end
        
        function train_by_ion(obj, trainVectors, trainData)
            warning('off', 'stats:classreg:learning:modifier:AdaBoostM1:modify:Terminate');
            
            obj.mClassifier = [];
            posSet = {'B', 'B', 'B', 'B', 'C', 'C', 'C', 'C'}; % 8
            negSet = {'C', 'Y', 'Z', 'O', 'B', 'Y', 'Z', 'O'};
            obj.mClassifier = [];
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
                    posSource = [posSource; [[trainData.(n).spectrumID]; [trainData.(n).peakID]]'];
                end
                
                for n = negIon
                    negX = [negX; trainVectors.(n)];
                    negSource = [negSource; [[trainData.(n).spectrumID]; [trainData.(n).peakID]]'];
                end
                
                X = [posX; negX];
                Y = [ones(size(posX,1), 1); zeros(size(negX,1), 1)];
                XInfo = [posSource; negSource];
                weights = [ones(size(posX,1),1)/size(posX,1); ones(size(negX,1),1)/size(negX,1)] / 2;
                
                temp = [trainVectors.massFeatures, trainVectors.massFeatures, trainVectors.massFeatures];
                obj.mClassifier.(p_v_n).X = X;
                obj.mClassifier.(p_v_n).Y = Y;
                obj.mClassifier.(p_v_n).XInfo = XInfo;
                obj.mClassifier.(p_v_n).BoostingNum = obj.mBoostNum;
                obj.mClassifier.(p_v_n).BoostingHoldout = obj.mHoldOut;
                
                % remove shared samples
                [~, ia] = unique( XInfo, 'rows' );
                
                if any(strcmp( p_v_n, use_simple_model ) )
                    obj.mClassifier.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
                    % obj.mClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'Bag', obj.mBoostNum, 'Tree', 'Holdout', obj.mHoldOut, 'Weights', weights(ia));
                    obj.mClassifier.(p_v_n).FeatureImportance = obj.mClassifier.(p_v_n).classifier.predictorImportance;
                else
                    obj.mClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'AdaBoostM1', obj.mBoostNum, 'Tree', 'Holdout', obj.mHoldOut, 'Weights', weights(ia));
                    %obj.mClassifier.(p_v_n).classifier = fitcensemble(X(ia,:), Y(ia,:), 'Method', 'AdaBoostM1', ...
                    %                                                'NumLearningCycles', obj.mBoostNum, 'Learners', 'tree', ...
                    %                                                'CrossVal', 'on', 'Holdout', obj.mHoldOut, 'Weights', weights(ia));
                    if isempty( obj.mClassifier.(p_v_n).classifier.Trained{1}.predictorImportance )
                        % obj.mClassifier.(p_v_n).classifier = fitctree(X(ia,:), Y(ia,:));
                        obj.mClassifier.(p_v_n).classifier = fitensemble(X(ia,:), Y(ia,:), 'Bag', obj.mBoostNum, 'Tree', ...
                            'Holdout', obj.mHoldOut, 'Weights', weights(ia), 'Type', 'classification');
                        %obj.mClassifier.(p_v_n).classifier = fitcensemble(X(ia,:), Y(ia,:), 'Method', 'Bag', ...
                        %                                            'NumLearningCycles', obj.mBoostNum, 'Learners', 'tree', ...
                        %                                            'CrossVal', 'on', 'Holdout', obj.mHoldOut, 'Weights', weights(ia));
                    end
                    obj.mClassifier.(p_v_n).FeatureImportance = obj.mClassifier.(p_v_n).classifier.Trained{1}.predictorImportance;
                    obj.mClassifier.(p_v_n).SelectedFeatureIdx = find( obj.mClassifier.(p_v_n).FeatureImportance(1:end-2) > 0 );
                    obj.mClassifier.(p_v_n).SelectedFeatureMass = temp( obj.mClassifier.(p_v_n).SelectedFeatureIdx );
                end
                disp( obj.mClassifier.(p_v_n) );
                % disp( num2cell([LOOResult.(p_v_n).SelectedFeatureIdx; LOOResult.(p_v_n).SelectedFeatureMass]) );
                
            end % round
            
            warning('on', 'stats:classreg:learning:modifier:AdaBoostM1:modify:Terminate');
        end
    end
    
    methods(Static)
        function [testVectors, dataSignals] = extract_data( spectra, massFeatures, massBound, massAccuracy, useOriginalPeak, includeREM )
            % Extract feature vectors of peaks with reconstruction results.
            if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 105]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            if nargin < 5 || isempty(useOriginalPeak), useOriginalPeak = 0; end % do not use the complementary peaks.
            if nargin < 6 || isempty(includeREM), includeREM = 0; end % include reducing end modification as a feature.
            
            allREM = {'', CMass.cReducingEndModification_O18, CMass.cReducingEndModification_Deuterium, ...
                CMass.cReducingEndModification_Reduced, CMass.cReducingEndModification_Aminopyridine, ...
                CMass.cReducingEndModification_PRAGS};
            
            dataSignals.B = []; % data of B-ions
            dataSignals.C = []; % data of C-ions
            
            for s = 1 : length(spectra)
                specU = spectra(s);
                if isempty( specU.mPeaks(end).mInferredSuperSet )
                    continue;
                end
                
                hasRequestREM = find( strcmp( specU.mReducingEndModification, allREM ) );
                
                peakMasses = [specU.mPeaks.mMass];
                peakComplements = [specU.mPeaks.mComplement] < 0;
                peakAvailable = zeros(1, length(specU.mPeaks));
                for p = 1 : length(specU.mPeaks)-1
                    if ~isempty(specU.mPeaks(p).mInferredFormulas)
                        peakAvailable(p) = 1;
                        
                        if useOriginalPeak && specU.mPeaks(p).mComplement < 0
                            peakAvailable(p) = 0;
                            peakAvailable(-specU.mPeaks(p).mComplement) = 1;
                        end
                    end
                end
                peakIntensities = CIonClassifier.standardize_intensity( specU.mPeaks, 'log', find(peakAvailable == 1) );
                
                for k = 1 : length(specU.mPeaks)
                    if isempty( specU.mPeaks(k).mInferredFormulas )
                        continue;
                    end
                    for TSS = specU.mPeaks(k).mInferredSuperSet
                        if isempty(TSS.mFormulas), continue; end % Modified: 2019/02/15
                        type = TSS.mTargetPeaks(2, ( TSS.mTargetPeaks(1,:) == k ));
                        switch type
                            case {1, 2}
                                pMass = peakMasses(k);
                            case {11, 21}
                                pMass = peakMasses(k) + CMass.H;
                            case {12, 22}
                                pMass = peakMasses(k) + CMass.H2;
                        end
                        switch type
                            case {1, 11, 12}
                                iontype = 'B';
                            case {2, 21, 22}
                                iontype = 'C';
                        end
                        temp = abs(peakMasses - pMass);
                        idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                        disp([iontype, ' ', num2str(k)]);
                        dataSignals.(iontype)(end+1).mass = pMass;
                        dataSignals.(iontype)(end).peakID = k;
                        dataSignals.(iontype)(end).complementPeakID = specU.mPeaks(k).mComplement;
                        dataSignals.(iontype)(end).intensity = peakIntensities(k);
                        dataSignals.(iontype)(end).context = [peakMasses(idx) - pMass; peakIntensities(idx); peakComplements(idx)];
                        dataSignals.(iontype)(end).spectrumID = s;
                        dataSignals.(iontype)(end).glycanMass = specU.mPrecursor;
                        dataSignals.(iontype)(end).REM = hasRequestREM;
                        dataSignals.(iontype)(end).spectrum = specU;
                    end
                end
            end
            
            disp( 'Converting into feature vectors ...' );
            testVectors.B = [];
            testVectors.C = [];
            
            numFeatures = length( massFeatures );
            for iontype = 'BC'
                num = length( dataSignals.(iontype) );
                testVectors.(iontype) = zeros(num, numFeatures * 2 + 3);
                % testVectors.(iontype) = zeros(num, numFeatures * 3 + 3); % complement flag. Not helpful.
                
                for k = 1 : num
                    testVectors.(iontype)(k, end) = dataSignals.(iontype)(k).REM;
                    testVectors.(iontype)(k, end-1) = dataSignals.(iontype)(k).glycanMass;
                    testVectors.(iontype)(k, end-2) = dataSignals.(iontype)(k).mass;
                    for m = 1 : size( dataSignals.(iontype)(k).context, 2 )
                        [d, idx] = min( abs( massFeatures - dataSignals.(iontype)(k).context(1, m) ) );
                        step = 0;
                        if d < massAccuracy + 0.001
                            step = 1;
                            testVectors.(iontype)(k, idx) = 1; % has the mass feature
                            testVectors.(iontype)(k, idx + numFeatures) = dataSignals.(iontype)(k).context(2, m); % intensity value
                            % testVectors.(iontype)(k, idx + 2*numFeatures) = dataSignals.(iontype)(k).context(3, m); % complement flag
                        end
                        if step == 0
                            % disp(iontype, num2str(k), num2str(m));
                        end
                    end
                end
                
                if ~includeREM
                    testVectors.(iontype) = testVectors.(iontype)(:, 1:end-1);
                end
            end
            testVectors.massFeatures = massFeatures;
            disp( 'CIonClassifier.extract_data() Done' );
        end
        
        function [dataVectors, dataSignals, massFeatures] = prepare_training_data( spectra, exprMethod, massBound, massAccuracy, useOriginalPeak, includeREM )
            % Each spectrum should provide its glycan.
            if nargin < 2, exprMethod = []; end % for example, exprMethod = 'EED'
            if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 105]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            if nargin < 5 || isempty(useOriginalPeak), useOriginalPeak = 0; end % useOriginalPeak = 1 do not use the complementary peaks.
            if nargin < 6 || isempty(includeREM), includeREM = 0; end % include reducing end modification as a feature.
            
            dataSignals.B = []; % data of B-ions
            dataSignals.C = []; % data of C-ions
            dataSignals.Y = []; % data of Y-ions
            dataSignals.Z = []; % data of Z-ions
            dataSignals.O = []; % data of O-ions
            
            glycanNames = {};
            for s = 1 : length(spectra)
                specU = spectra(s);
                
                % filter by experiment method
                if ~isempty(exprMethod) && strcmp(specU.mExperimentMethod, exprMethod ) == 0
                    continue;
                end
                
                numbers = zeros(1, 4);
                peakMasses = [specU.mPeaks.mMass];
                peakComplements = [specU.mPeaks.mComplement] < 0;
                peakAvailable = zeros(1, length(specU.mPeaks));
                for p = 1 : length(specU.mPeaks)-1
                    if ~isempty(specU.mPeaks(p).mInferredFormulas)
                        peakAvailable(p) = 1;
                        
                        if useOriginalPeak && specU.mPeaks(p).mComplement < 0
                            peakAvailable(p) = 0;
                            peakAvailable(-specU.mPeaks(p).mComplement) = 1;
                        end
                    end
                end
                peakIntensities = CIonClassifier.standardize_intensity( specU.mPeaks, 'log', find(peakAvailable == 1) );
                % Features will be built from peakMasses, peakComplements, and peakIntensities.
                
                g = CGlycan( specU.mPermethylated, specU.mReducingEndModification );
                idx = strfind( specU.comment, '.' );
                name = strtrim( specU.comment(idx+1:end) );
                try
                    g.parse( name );
                catch e
                    disp( e ); disp( s );
                    continue;
                end
                glycanNames{end+1} = name;
                ions = [];
                [ions.Y, ions.B, ~, ~, ions.Z, ions.C] = g.cleave( 'Proton' ); % [Ys, Bs, Xs, As, Zs, Cs] = cleave( obj, specU.mMetal );
                allREM = {'', CMass.cReducingEndModification_O18, CMass.cReducingEndModification_Deuterium, ...
                    CMass.cReducingEndModification_Reduced, CMass.cReducingEndModification_Aminopyridine, ...
                    CMass.cReducingEndModification_PRAGS};
                hasRequestREM = find( strcmp( specU.mReducingEndModification, allREM ) );
                
                ionprops = [];
                ionTypes = {'B', 'C', 'Y', 'Z'};
                for it = 1 : 4
                    iontype = ionTypes{it};
                    
                    % consolidate masses to deal with minor precision problem.
                    ionprops.(iontype).masses = [ions.(iontype).mMass];
                    [uniqueMasses, ~, uID] = unique( floor( ionprops.(iontype).masses * 1000 ) / 1000 );
                    numUniqueMasses = length(uniqueMasses);
                    for kk = 1 : numUniqueMasses
                        idx = abs(ionprops.(iontype).masses - uniqueMasses(kk)) < massAccuracy;
                        uniqueMasses(kk) = mean(ionprops.(iontype).masses(idx));
                    end
                    ionprops.(iontype).uniqueMasses = uniqueMasses;
                    
                    % collect signals
                    ionprops.(iontype).uniqueMassAvailable = ones(1, numUniqueMasses);
                    for k = 1 : numUniqueMasses
                        ionM = uniqueMasses(k);
                        [minDiff, minIdx] = min( abs( peakMasses - ionM ) );
                        if minDiff < massAccuracy
                            peakAvailable( minIdx ) = 0;
                            ionprops.(iontype).uniqueMassAvailable(k) = 0;
                            temp = abs(peakMasses - ionM);
                            idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), minIdx );
                            dataSignals.(iontype)(end+1).mass = peakMasses(minIdx);
                            dataSignals.(iontype)(end).peakID = minIdx;
                            dataSignals.(iontype)(end).complementPeakID = specU.mPeaks(minIdx).mComplement;
                            dataSignals.(iontype)(end).intensity = peakIntensities(minIdx);
                            dataSignals.(iontype)(end).context = [peakMasses(idx)-ionM; peakIntensities(idx); peakComplements(idx)];
                            dataSignals.(iontype)(end).glycanMass = g.mMass;
                            dataSignals.(iontype)(end).spectrumID = s;
                            dataSignals.(iontype)(end).REM = hasRequestREM;
                            dataSignals.(iontype)(end).spectrum = spectra(s);
                            dataSignals.(iontype)(end).ions = ions.(iontype)(uID == k);
                            numbers(it) = numbers(it) + 1;
                        end
                    end
                end
                
                for k = 1 : length(specU.mPeaks)-1
                    if ~peakAvailable(k)
                        continue;
                    end
                    for TSS = specU.mPeaks(k).mInferredSuperSet
                        iontype = '';
                        type = TSS.mTargetPeaks(2, ( TSS.mTargetPeaks(1,:) == k ));
                        switch type %here: need +CMASS?
                            case {1, 2}
                                pMass = peakMasses(k);
                            case {11, 21}
                                pMass = peakMasses(k) + CMass.H;
                            case {12, 22}
                                pMass = peakMasses(k) + CMass.H2;
                        end
                        switch type
                            case {1, 11, 12}
                                iontype = 'B';
                            case {2, 21, 22}
                                iontype = 'C';
                        end
                        [minDiff, minIdx] = min( abs( pMass - ionprops.(iontype).uniqueMasses ) );
                        if minDiff < massAccuracy && ionprops.(iontype).uniqueMassAvailable(minIdx) % Use uniqueMassAvailable to avoid duplication
                            peakAvailable( minIdx ) = 0;
                            ionprops.(iontype).uniqueMassAvailable(minIdx) = 0;
                            temp = abs(peakMasses - pMass);
                            %2020/10/17 here
                            idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                            dataSignals.(iontype)(end+1).mass = pMass;
                            dataSignals.(iontype)(end).peakID = k;
                            dataSignals.(iontype)(end).complementPeakID = specU.mPeaks(k).mComplement;
                            dataSignals.(iontype)(end).intensity = peakIntensities(k);
                            dataSignals.(iontype)(end).context = [peakMasses(idx) - pMass; peakIntensities(idx); peakComplements(idx)];
                            dataSignals.(iontype)(end).spectrumID = s;
                            dataSignals.(iontype)(end).glycanMass = g.mMass;
                            dataSignals.(iontype)(end).REM = hasRequestREM;
                            dataSignals.(iontype)(end).spectrum = spectra(s);
                            dataSignals.(iontype)(end).ions = ions.(iontype)(uID == minIdx);
                        end
                    end
                end
                
                disp( [name, ': ', sprintf( '%d ', numbers )] );
                if numbers(1) == 0, disp( ['Missing B: ', name]); end
                if numbers(2) == 0, disp( ['Missing C: ', name]); end
                if numbers(3) == 0, disp( ['Missing Y: ', name]); end
                if numbers(4) == 0, disp( ['Missing Z: ', name]); end
                
                for k = 1 : length(specU.mPeaks)-1
                    if peakAvailable(k)
                        peakAvailable(k) = 0;
                        temp = abs( peakMasses - peakMasses(k) );
                        idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                        dataSignals.O(end+1).mass = peakMasses(k);
                        dataSignals.O(end).peakID = k;
                        dataSignals.O(end).complementPeakID = specU.mPeaks(k).mComplement;
                        dataSignals.O(end).intensity = peakIntensities(k);
                        dataSignals.O(end).context = [peakMasses(idx) - peakMasses(k); peakIntensities(idx); peakComplements(idx)];
                        dataSignals.O(end).spectrumID = s;
                        dataSignals.O(end).glycanMass = g.mMass;
                        dataSignals.O(end).REM = hasRequestREM;
                        dataSignals.O(end).spectrum = spectra(s);
                    end
                end
            end
            
            disp( 'Converting into feature vectors ...' );
            dataVectors.B = [];
            dataVectors.C = [];
            dataVectors.Y = [];
            dataVectors.Z = [];
            dataVectors.O = [];
            
            % build feature set
            masses = [];
            for iontype = 'BCYZO'
                for k = 1 : length( dataSignals.(iontype) )
                    masses = [masses, dataSignals.(iontype)(k).context(1,:)];
                end
            end
            masses = round( masses * 1000 ) / 1000;
            masses = unique(masses);
            massFeatures = merge_masses(masses', massAccuracy, 1);
            numFeatures = length( massFeatures );
            
            for iontype = 'BCYZO'
                num = length( dataSignals.(iontype) );
                dataVectors.(iontype) = zeros(num, numFeatures * 2 + 3);
                % dataVectors.(iontype) = zeros(num, numFeatures * 3 + 3); % complementary flag, not helpful
                
                for k = 1 : num
                    dataVectors.(iontype)(k, end) = dataSignals.(iontype)(k).REM;
                    dataVectors.(iontype)(k, end-1) = dataSignals.(iontype)(k).glycanMass;
                    dataVectors.(iontype)(k, end-2) = dataSignals.(iontype)(k).mass;
                    for m = 1 : size( dataSignals.(iontype)(k).context, 2 )
                        [d, idx] = min( abs( massFeatures - dataSignals.(iontype)(k).context(1, m) ) );
                        if d < massAccuracy + 0.001
                            dataVectors.(iontype)(k, idx) = 1; % has the mass feature
                            dataVectors.(iontype)(k, idx + numFeatures) = dataSignals.(iontype)(k).context(2, m); % intensity value
                            % dataVectors.(iontype)(k, idx + numFeatures*2) = dataSignals.(iontype)(k).context(3, m); % complementary flag
                        else
                            disp( ['Missing mass features! ', iontype, ': ', num2str(k)] );
                        end
                    end
                end
                
                if ~includeREM
                    dataVectors.(iontype) = dataVectors.(iontype)(:, 1:end-1);
                end
            end
            dataVectors.massFeatures = massFeatures;
            disp( 'CIonClassifier.prepare_training_data() Done' );
        end
        
        function [filteredSignals, filteredVectors] = filter_BC_by_rootmono( dataSignals, dataVectors, rootmonos )
        % [filteredSignals, filteredVectors] = filter_BC_by_rootmono( dataSignals, dataVectors, rootmonos )
            for iontype = 'BC'
                N = length(dataSignals.(iontype));
                flag = zeros(1, N);
                for k = 1 : N
                    for m = 1 : length( dataSignals.(iontype)(k).ions )
                        if any(strcmp(dataSignals.(iontype)(k).ions(m).mStem(1).mSymbol, rootmonos))
                            flag(k) = 1;
                            break;
                        end
                    end
                end
                flag = flag > 0;
                filteredVectors.(iontype) = dataVectors.(iontype)(flag,:);
                filteredSignals.(iontype) = dataSignals.(iontype)(flag);
            end
            disp( 'CIonClassifier.filter_data_by_rootmono() Done' );
        end
        
        function [zscores, intensities, m, s] = standardize_intensity( peaks, transform, peakIDs )
            if nargin < 2, transform = ''; end
            if nargin < 3, peakIDs = []; end
            
            intensities = [peaks.mIntensity];
            switch transform
                case 'log'
                    intensities = log(intensities);
                case 'sqrt'
                    intensities = sqrt(intensities);
            end
            
            % [m, std] = robust_mean_std( intensities( [peaks.mComplement] == 0 ) );
            if isempty( peakIDs )
                [m, s] = robust_mean_std( intensities(1:end-1) );
                zscores = (intensities - m)/s;
            else
                [m, s] = robust_mean_std( intensities(peakIDs) );
                zscores = (intensities - m)/s;
            end
        end
        
    end % methods(Static)
    
end