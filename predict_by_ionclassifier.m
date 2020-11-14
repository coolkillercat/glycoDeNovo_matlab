function [results, ICScores] = predict_by_ionclassifier( spectra, dataPeaks, dataVectors, icPool )
% This function is written for LOO setting. 

num = length( spectra );
ICScores = cell(1, num);
results = {'SpecID', 'Name', 'Formula', 'Fragmentation', 'Precursor', 'ReducingEndModification', ...
           'Metal', '# Peaks', '# Interpretable Peaks', '# Reconstructed Peaks', '# Candidates', ...
           'Rank by Supporting Peak Number', 'Rank by IonClassifier', 'Best IC', '2nd IC', 'File'};
scoreMapBIon = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
scoreMapCIon = containers.Map('KeyType', 'uint32', 'ValueType', 'double');
spectraMap = {};
for ion = 'BCYZO'
    for k = 1 : length(dataPeaks.(ion))
        spectraMap{dataPeaks.(ion)(k).glycanID} = dataPeaks.(ion)(k).spectrum;
        key = uint32(dataPeaks.(ion)(k).glycanID * 10000 + dataPeaks.(ion)(k).peakID);
        
        % Calculate its score of being a B-ion
        ws = zeros(1, 4);
        [~, b] = icPool.B_v_C.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.B_v_C.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.B_v_C.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.B_v_Y.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.B_v_Y.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.B_v_Y.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.B_v_Z.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.B_v_Z.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.B_v_Z.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.B_v_O.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.B_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.B_v_O.classifier.NumTrainedPerFold;
        end
        scoreMapBIon(key) = min(ws);

        % Calculate its score of being a C-ion
        [~, b] = icPool.C_v_B.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.C_v_B.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.C_v_B.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.C_v_Y.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.C_v_Y.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.C_v_Y.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.C_v_Z.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.C_v_Z.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.C_v_Z.classifier.NumTrainedPerFold;
        end
        [~, b] = icPool.C_v_O.classifier.Trained{1}.predict( dataVectors.(ion)(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.C_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.C_v_O.classifier.NumTrainedPerFold;
        end
        scoreMapCIon(key) = min(ws);
    end
end

for s = 1 : num
    if isempty(spectra(s).mPeaks(end).mInferredSuperSet)
        continue;
    end

    glycanID = -1;
    for g = 1 : length(spectraMap)
        if spectra(s) == spectraMap{g}
            glycanID = g;
        end
    end
    assert( glycanID > 0 );
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
            
            key = uint32(glycanID * 10000 + peakID);
            
            if type == 1 || type == -1
                weights(m) = scoreMapBIon(key);
            elseif type == 2 || type == -2
                weights(m) = scoreMapCIon(key);
            end
            ICScores{s}{t}(4, m) = weights(m);
        end        
        tp.mScore = sum(weights);
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
    name = [];
    if ~isempty(idx)
        name = spectra(s).comment(1:idx-1);
        formula = spectra(s).comment(idx+2:end);
    else
        formula = spectra(s).comment;
    end
    formula = strtrim( formula );
    formula = strrep( formula, '][', '] [' );
    
    if isempty( name )
        name = spectra(s).filename;
        idx = strfind( name, '\' );
        name = spectra(s).filename(idx(end)+1:end);
        idx = strfind( name, '.' );
        if length( idx ) > 1
            name = name(idx(1)+1:idx(2)-1);
        else
            name = strrep( name(5:end), '.mat', '' );
        end
    end
    
    % Get groundtruth
    g = CGlycan;
    g.parse( formula );
    g.order_branch_by_name();
    
    % get rank information
    rankInfo = spectra(s).get_rank( g.mInferredFormula );
    
    if rankInfo.rankBySupportingPeak == 0, continue; end
    
    numPeaks = [num2str(length(spectra(s).mPeaks)), '(', num2str(sum( [spectra(s).mPeaks.mComplement] < 0)), ')'];
    byPeaks = [num2str(rankInfo.rankBySupportingPeak), '(', num2str(rankInfo.corankBySupportingPeak-1), ')'];
    byIonClassifier = [num2str(rankInfo.rankByIonclassifier), '(', num2str(rankInfo.corankByIonclassifier-1), ')'];
    disp( [byPeaks, ', ', byIonClassifier, ', ', formula, ': ', spectra(s).filename] );
    results(end+1,:) = {s, name, formula, spectra(s).mExperimentMethod, spectra(s).mPrecursor, ...
        spectra(s).mReducingEndModification, spectra(s).mMetal, ...
        numPeaks, numInterpretable, numReconstructed, ...
        length(spectra(s).mPeaks(end).mInferredFormulas), ...
        byPeaks, byIonClassifier, rankInfo.topIonCalssifierScores(end), ...
        rankInfo.topIonCalssifierScores(end-1), spectra(s).filename };
end
%disp( 'Finished' );
