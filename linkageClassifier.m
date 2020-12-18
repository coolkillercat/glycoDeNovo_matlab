classdef linkageClassifier
    
    properties
        mMassFeatures = [];
        mMassAccuracy = 0.005;
        mClassifier = [];
        mBoostNum = 20;
        mClassifierMap; %one Map for one ion
        mMethods = ["SVM"];
        linkageList = ["HexHex", "HexNAcHex", "HexHexNAc"];
    end
    
    methods
        %TODO: estimate linkage using spectrumdata; generate labels using
        %annotation data; calculate loss (ask Zizhang about how to calculate p value?)
        
        function classifier = train(~, X, Y, method, hp)
            classifier = trainMethods.train_by_method(X, Y, method, hp);
        end
        
        function obj = train_by_ion(obj, dataVectors, dataSources)
            obj.mClassifierMap = containers.Map;
            vs = nchoosek(keys(dataVectors),2);
            for m = obj.mMethods %Edit method here
                for v = 1 : length(vs)
                    pos = vs{v,1};
                    neg = vs{v,2};
                    posX = dataVectors(pos);
                    posSource = dataSources(pos);
                    negX = dataVectors(neg);
                    negSource = dataSources(neg);
                    X = [posX; negX];
                    Y = [ones(size(posX,1), 1); zeros(size(negX,1), 1)];
                    XInfo = [posSource, negSource];
                    weights = [ones(size(posX,1),1)/size(posX,1); ones(size(negX,1),1)/size(negX,1)] / 2;
                    hp = cell(2,1);
                    hp{1} = obj.mBoostNum;
                    hp{2} = weights;
                    currclassifier = obj.train(X, Y, m, hp);
                    trainLabel = predict(currclassifier.classifier, X);
                    if strcmp(m, 'randomForest') == 1
                         trainLabel = str2double(trainLabel);
                    end
                    trainError = sum(abs(trainLabel - Y)) / length(Y);
                    currclassifier.trainError = trainError;
                    currclassifier.XInfo = XInfo;
                    obj.mClassifierMap([pos, neg]) = currclassifier;
                end
            end
        end
        
        %dataVectors should have 1 ion
        function linkageScoreMap = predict_linkages(obj, dataVectors, dataSources, linkageList)
            if nargin < 4, linkageList = obj.linkageList; end
            linkageScoreMap = containers.Map;
            %========================================
            avlClassifiers = zeros(length(linkageList), length(linkageList));
            for linkage = 1:length(linkageList)
                for vsLinkage = 1:length(linkageList)
                    key = strcat(linkageList(linkage), linkageList(vsLinkage));
                    if isKey(obj.mClassifierMap, key)
                        linkageScoreMap(linkageList(linkage)) = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
                        linkageScoreMap(linkageList(vsLinkage)) = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
                        avlClassifiers(linkage, vsLinkage) = 1;
                        avlClassifiers(vsLinkage, linkage) = -1;
                    end
                end
            end
            %==========================================
            for k = 1 : size(dataVectors,1)
                key = dataSources(k).spectrumID * 10000 + dataSources(k).peakID;
                for linkage = 1 : length(linkageList)
                    if ~isKey(linkageScoreMap, linkageList(linkage))
                        continue;
                    end
                    vsLinkages = find(avlClassifiers(linkage,:)~=0);
                    ws = zeros(1, length(vsLinkages));
                    i = 0;
                    for vsLinkage = vsLinkages
                        i = i + 1;
                        if avlClassifiers(linkage, vsLinkage) == 1
                            currClassifier = obj.mClassifierMap(strcat(linkageList(linkage), linkageList(vsLinkage)));
                            [tttt, b] = currClassifier.classifier.predict( dataVectors(k,:) );
                            ws(i) = b(2);
                        else
                            currClassifier = obj.mClassifierMap(strcat(linkageList(vsLinkage), linkageList(linkage)));
                            [tttt, b] = currClassifier.classifier.predict( dataVectors(k,:) );
                            ws(i) = b(1);
                        end
                    end
                    map = linkageScoreMap(linkageList(linkage));
                    map(key) = min(ws);%?Dont know if it works
                end
            end
        end
        
        function linkageScoreMap = rank_candidates(obj, specSet)
            [dataVectors, dataSources] = obj.extract_data(specSet);
            linkageScoreMap = obj.predict_linkages(dataVectors.B, dataSources.B);%TODO
        end
        
        function  [accurateNum, accuracy, wrongList, wrongData] = testScoreMap(obj, linkageScoreMap, specSetTrain, specSetTest)
            linkages = keys(linkageScoreMap);
            peakIDs = keys(linkageScoreMap(linkages{1}));
            accurateNum = 0;
            wrongList = [];
            wrongData = [];
            for p = 1 : length(peakIDs)
                peakID = peakIDs{p};
                specTrain = specSetTrain{floor(peakID/10000)};
                specTest = specSetTest{floor(peakID/10000)};
                peakTrain = specTrain.mPeaks(mod(peakID, 10000));
                [~, minIdx] = min(abs([specTest.mPeaks.mRawMZ] - peakTrain.mMass));
                peakTest = specTest.mPeaks(minIdx);
                maxScore = 0;
                maxLinkage = [];
                scores = [];
                for l = 1:length(linkages)
                    linkage = linkages{l};
                    lmap = linkageScoreMap(linkage);
                    score = lmap(peakID);
                    scores = [scores, score];
                    if score > maxScore
                        maxScore = score;
                        maxLinkage = linkage;
                    end
                end
                [RE, NRE] = linkageClassifier.split_linkage(maxLinkage);
                peakTrain.RE = RE;
                peakTrain.NRE = NRE;
                if strcmp(maxLinkage, [peakTest.RE, peakTest.NRE]) == 1
                    accurateNum = accurateNum + 1;
                else
                    if strcmp(peakTest.RE, '') == 0 || 1
                        wrongList = [wrongList, peakTrain];
                        wrongData = [wrongData; [double(peakID), scores]];
                    end
                end
            end
            wrongData = [["peakIDs", linkages]; wrongData];
            accuracy = accurateNum / length(peakIDs);
        end
        
        function massFeatures = calculate_mass_feature(obj, specSet, masses, massBound, massAccuracy)
            if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 30]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            h1 = waitbar(0, 'spectrum');
            h2 = waitbar(0, 'peak');
            massFeatures = merge_masses(masses', massAccuracy, 1);
            for s = 1 : length(specSet)
                spectra = specSet{s};
                waitbar(s/length(specSet), h1 , ['spectrum progress...', num2str(100*s/length(specSet)), '%']);
                for k = 1 : length(spectra.mPeaks)
                    if (mod(k + 1, 10) == 0) 
                        waitbar(k/length(spectra.mPeaks), h2, ['peak progress...', num2str(100*k/length(spectra.mPeaks)), '%'])
                    end
                    peak = spectra.mPeaks(k);
                    masslen = length(massFeatures);
                    peak.mLinkageVector = [peak.mRawMZ, zeros(1, 2 * masslen)];
                    featureAvailable = ones(1, length(massFeatures));
                    neighbour = spectra.mPeaks(k).mNeighbour;
                    for kk = 1 : length(neighbour(1,:))
                        [d, idx] = min( abs( massFeatures - neighbour(1,kk) ) );
                        if d < massAccuracy + 0.001 && featureAvailable(idx)
                            peak.mLinkageVector(idx + 1) = 1;
                            peak.mLinkageVector(idx + 1 + masslen) = spectra.mPeaks(neighbour(2,kk)).mZscore;
                            spectra.mPeaks(k).mFeature = [spectra.mPeaks(k).mFeature, massFeatures(idx)];
                            spectra.mPeaks(k).mFeatureIntensity = [spectra.mPeaks(k).mFeatureIntensity, spectra.mPeaks(neighbour(2,kk)).mZscore];
                        end
                    end
                end
            end
            delete(h1);
            delete(h2);
        end
        
        function [dataVectors, dataSources] = extract_data(obj, specSet, massBound, massAccuracy)
            if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 30]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            peakAvailable = cell(length(specSet), 1);
            for s = 1: length(specSet)
                spectra = specSet{s};
%                 if ~spectra.mProtonated %test: always protonate
%                     spectra.protonate();
%                     spectra.add_complementary_ions();
%                     spectra.merge_peaks( 0.0025 );
% %                     m2c = CMass2Composition;
% %                     m2c.load( 'm2c.mat' );
% %                     m2c.mCheckMinus2H = 1;
% %                     if spectra.mMassAccuracy > 0
% %                         m2c.mMassAccuracyPPM = spectra.mMassAccuracy;
% %                     else
% %                         m2c.mMassAccuracyPPM = 5;
% %                     end
% %                     m2c.set_reducing_end_modification( spectra.mReducingEndModification );
% %                     m2c.set_permethylation( spectra.mPermethylated );
% %                     hypoSpecs = m2c.correct_spectrum( spectra );
%                     reconstructor = CGlycoDeNovo( 5, []); % 5ppm accuracy
%                     reconstructor.mCheckMinusH = 0;
%                     reconstructor.interpret_peaks( spectra );
%                     reconstructor.reconstruct_formulas();
%                 end
                context = {};
                peakMasses = [spectra.mPeaks.mMass];
                neighbour = cell(length(peakMasses), 1);
                peakAvailable{s} = ones(1, length(peakMasses)); %should be zeros
                for p = 1 : length(spectra.mPeaks)-1
                    if ~isempty(spectra.mPeaks(p).mInferredFormulas)  %?:no peaks available
                        peakAvailable{s}(p) = 1;
                    end
                end
                for k = 1 : length(spectra.mPeaks)
                    if ~peakAvailable{s}(k), continue; end
                    pMass = peakMasses(k);
                    temp = abs(peakMasses - pMass);
                    idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                    neighbour = [peakMasses(idx) - pMass; idx];
                    spectra.mPeaks(k).mNeighbour = neighbour;
                end
                peakIntensities = linkageClassifier.standardize_intensity( spectra.mPeaks, 'log', find(peakAvailable{s} == 1) );
                
            end
            h1 = waitbar(0, 'spectrum');
            h2 = waitbar(0, 'peak');
            massFeatures = obj.mMassFeatures;
            for s = 1 : length(specSet)
                spectra = specSet{s};
                waitbar(s/length(specSet), h1 , ['spectrum progress...', num2str(100*s/length(specSet)), '%']);
                for k = 1 : length(spectra.mPeaks)
                    if ~peakAvailable{s}(k), continue; end
                    
                    if (mod(k + 1, 10) == 0) 
                        waitbar(k/length(spectra.mPeaks), h2, ['peak progress...', num2str(100*k/length(spectra.mPeaks)), '%'])
                    end
                    peak = spectra.mPeaks(k);
                    masslen = length(massFeatures);
                    peak.mLinkageVector = [peak.mMass, zeros(1, 2 * masslen)];
                    featureAvailable = ones(1, length(massFeatures));
                    neighbour = spectra.mPeaks(k).mNeighbour;
                    for kk = 1 : length(neighbour(1,:))
                        [d, idx] = min( abs( massFeatures - neighbour(1,kk) ) );
                        if d < massAccuracy + 0.001 && featureAvailable(idx)
                            peak.mLinkageVector(idx + 1) = 1;
                            peak.mLinkageVector(idx + 1 + masslen) = spectra.mPeaks(neighbour(2,kk)).mZscore;
                            spectra.mPeaks(k).mFeature = [spectra.mPeaks(k).mFeature, massFeatures(idx)];
                            spectra.mPeaks(k).mFeatureIntensity = [spectra.mPeaks(k).mFeatureIntensity, spectra.mPeaks(neighbour(2,kk)).mZscore];
                        end
                    end
                end
            end
            delete(h1);
            delete(h2);
            dataVectors.B = [];
            dataVectors.C = [];
            dataSources.B = [];
            dataSources.C = [];
            for s = 1: length(specSet)
                spec = specSet{s};
                for p = 1 : length(spec.mPeaks)
                    if ~peakAvailable{s}(p), continue; end
                    peak = spec.mPeaks(p);
                    for TSS = peak.mInferredSuperSet
                       if isempty(TSS.mFormulas), continue; end
                        type = TSS.mTargetPeaks(2, ( TSS.mTargetPeaks(1,:) == p ));
                        switch type
                            case {1, 2}
                                pMass = peakMasses(p);
                            case {11, 21}
                                pMass = peakMasses(p) + CMass.H;
                            case {12, 22}
                                pMass = peakMasses(p) + CMass.H2;
                        end
                        switch type
                            case {1, 11, 12}
                                iontype = 'B';
                            case {2, 21, 22}
                                iontype = 'C';
                        end
                        dataSources.(iontype)(end+1).peak = peak;
                        dataSources.(iontype)(end).context = peak.mNeighbour;
                        dataSources.(iontype)(end).peakID = p;
                        dataSources.(iontype)(end).spectrumID = s;
                        dataVectors.(iontype) = [dataVectors.(iontype);peak.mLinkageVector];
                    end
                end
            end
        end
        
        
    end
    
    methods(Static)
        
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
            for p = 1 : length(peaks)
                    peaks(p).mZscore = zscores(p);
            end
        end
        
        function [dataVectors, dataSources] = prepare_training_data(traindata)
            dataVectors = containers.Map;
            dataSources = containers.Map;
            for k = 1 : length(traindata)
                peak = traindata(k);
                linkage = [peak.RE, peak.NRE];
                if ~isKey(dataVectors, linkage)
                    dataVectors(linkage) = [peak.mLinkageVector];
                    dataSources(linkage) = [k];
                else
                    dataVectors(linkage) = [dataVectors(linkage); peak.mLinkageVector];
                    dataSources(linkage) = [dataSources(linkage), k];
                end
            end
        end
        
        function masses = calculate_masses(spectra,orimasses, massBound, massAccuracy)
             if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 30]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            context = {};
            peakMasses = [spectra.mPeaks.mRawMZ];
            neighbour = cell(length(peakMasses), 1);
            peakAvailable = ones(1, length(peakMasses));
            [uniqueMasses, ~, uID] = unique( floor( peakMasses * 1000 ) / 1000 );
            numUniqueMasses = length(uniqueMasses);
            for kk = 1 : numUniqueMasses
                  idx = abs(peakMasses - uniqueMasses(kk)) < massAccuracy;
                  uniqueMasses(kk) = mean(peakMasses(idx));
            end
            uniqueMassAvailable = ones(1, numUniqueMasses);
            for k = 1 : numUniqueMasses
            	ionM = uniqueMasses(k);
            	[minDiff, minIdx] = min( abs( peakMasses - ionM ) );
                	if minDiff < massAccuracy
                        peakAvailable( minIdx ) = 0;
                    	uniqueMassAvailable(k) = 0;
                    	temp = abs(peakMasses - ionM);
                    	idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), minIdx );
                    	context{end+1} = [peakMasses(idx)-ionM; idx];
               	end
            end
            for k = 1 : length(peakMasses)
                pMass = peakMasses(k);
                [minDiff, minIdx] = min( abs( pMass - uniqueMasses ) );
                temp = abs(peakMasses - pMass);
                idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                neighbour{k} = [peakMasses(idx) - pMass; idx];
                spectra.mPeaks(k).mNeighbour = neighbour{k};
                if minDiff < massAccuracy && uniqueMassAvailable(minIdx) % Use uniqueMassAvailable to avoid duplication
                	peakAvailable( minIdx ) = 0;
                	uniqueMassAvailable(minIdx) = 0;
                    context{end+1} = neighbour{k};
               end
            end
            masses = orimasses;
            for k = 1 : length(context)
                masses = [masses, context{k}(1,:)];
            end
            masses = round( masses * 1000 ) / 1000;
            masses = unique(masses);
        end
        
        function getLinkageScore(spec, LSM)
            if isempty(LSM.B), LSM.B = containers.Map(); end
            if isempty(LSM.C), LSM.C = containers.Map(); end
            for p = 1:length(spec.mPeaks)
                peak = spec.mPeaks(p);
                if isempty(peak.mInferredSuperSet), continue; end
                if linkageClassifier.isMono(peak)
                    peak.mInferredSuperSet.mTopologies.mLinkageScore = 0;
                    peak.mInferredSuperSet.mLinkageScore = 0;
                    peak.mInferredSuperSet.mTopologySets.mLinkageScore = 0;
                end
                for TSS = peak.mInferredSuperSet
                    for TS = TSS.mTopologySets
                        RE = TS.mRootMono.mClass;
                        for STSS = TS.Sources
                            cache = [];
                            if isempty(STSS), continue; end
                            for STS = STSS.mTopologySets
                                NRE = STS.mRootMono.mClass;
                                SPeak = STS.mTopologies(1).mSupportPeaks(end);
                                iontype = STS.mType;
                                imap = LSM.(iontype);
                                linkage = [RE, NRE];
                                lmap = imap(linkage);
                                addScore = lmap(SPeak);
                                cache = [cache, addScore + STS.mLinkageScore];
                            end
                        end
                    end
                end
            end
        end
        
        function linkageList = mono2linkage(monoList)
            vs = nchoosek(monoList,2);
            linkageList = [];
            for v = 1:size(vs,1)
                linkageList = [linkageList, strcat(vs(v,1),vs(v,2)) , strcat(vs(v,2), vs(v,1))];
            end
            for v = 1:length(monoList)
                linkageList = [linkageList, strcat(monoList(v), monoList(v))];
            end
        end
        
        function [RE, NRE] = split_linkage(linkage)
            ch = linkage(1);
            RL = 3;
            if strcmp(ch, 'H') == 1
                if strcmp(ch, 'HexNAc') == 1
                    RL = 6;
                elseif strcmp(ch, 'HexA') == 1
                        RL = 4;
                end
            elseif strcmp(ch, 'N') == 1
                RL = 6;
            end
            RE = linkage(1:RL);
            NRE = linkage((RL+1):end);
        end
        
        function label = isMono(peak)
            if isempty(peak.mInferredSuperSet), label = 0; return; end
            label = 1;
            for TSS = peak.mInferredSuperSet
                if length(TSS.mTopologies) == 1 && sum(TSS.mTopologies.mCompositionCount) == 1
                else 
                    label = 0;
                end
            end
        end
        
    end
    
end

