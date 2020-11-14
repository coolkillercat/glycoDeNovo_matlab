function [results, ICScores] = predict_by_ionclassifier( spectra, dataPeaks, dataVectors, icPool )

num = length( spectra );
ICScores = cell(1, num);
results = {'SpecID', 'Name', 'Formula', 'Fragmentation', 'Precursor', 'ReducingEndModification', ...
           'Metal', '# Peaks', '# Interpretable Peaks', '# Reconstructed Peaks', '# Candidates', ...
           'Rank by Supporting Peak Number', 'Rank by IonClassifier', 'File'};
for s = 1 : num
    if isempty(spectra(s).mPeaks(end).mInferredSuperSet)
        continue;
    end
    if strcmp(spectra(s).mExperimentMethod, 'EED' ) == 0
        continue;
    end

    scores.B = zeros(size( dataVectors.B, 1 ),4);
    for k = 1 : size( dataVectors.B, 1 )
        ws = zeros(1, 4);

        [a, b] = icPool.B_v_C.classifier.Trained{1}.predict( dataVectors.B(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.B_v_C.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.B_v_C.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.B_v_Y.classifier.Trained{1}.predict( dataVectors.B(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.B_v_Y.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.B_v_Y.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.B_v_Z.classifier.Trained{1}.predict( dataVectors.B(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.B_v_Z.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.B_v_Z.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.B_v_O.classifier.Trained{1}.predict( dataVectors.B(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.B_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.B_v_O.classifier.NumTrainedPerFold;
        end
        
        %scores.B(k) = min(ws);
        scores.B(k,:) = ws;
    end

    scores.C = zeros(size( dataVectors.C, 1 ),4);
    for k = 1 : size( dataVectors.C, 1 )
        ws = zeros(1, 4);

        [a, b] = icPool.C_v_B.classifier.Trained{1}.predict( dataVectors.C(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.C_v_B.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.C_v_B.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.C_v_Y.classifier.Trained{1}.predict( dataVectors.C(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.C_v_Y.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.C_v_Y.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.C_v_Z.classifier.Trained{1}.predict( dataVectors.C(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.C_v_Z.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.C_v_Z.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.C_v_O.classifier.Trained{1}.predict( dataVectors.C(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.C_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.C_v_O.classifier.NumTrainedPerFold;
        end
        
        %scores.C(k) = min(ws);
        scores.C(k,:) = ws;
    end
    
    scores.Y = zeros(size( dataVectors.Y, 1 ), 4);
    for k = 1 : size( dataVectors.Y, 1 )
        ws = zeros(1, 4);

        [a, b] = icPool.Y_v_B.classifier.Trained{1}.predict( dataVectors.Y(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.Y_v_B.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.Y_v_B.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Y_v_C.classifier.Trained{1}.predict( dataVectors.Y(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.Y_v_C.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.Y_v_C.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Y_v_Z.classifier.Trained{1}.predict( dataVectors.Y(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.Y_v_Z.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.Y_v_Z.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Y_v_O.classifier.Trained{1}.predict( dataVectors.Y(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.Y_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.Y_v_O.classifier.NumTrainedPerFold;
        end
        
        %scores.Y(k) = min(ws);
        scores.Y(k,:) = ws;
    end
    
    scores.Z = zeros(size( dataVectors.Z, 1 ), 4);
    for k = 1 : size( dataVectors.Z, 1 )
        ws = zeros(1, 4);

        [a, b] = icPool.Z_v_B.classifier.Trained{1}.predict( dataVectors.Z(k,:) );
        ws(1) = b(2);
        if strcmp( icPool.Z_v_B.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.Z_v_B.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Z_v_C.classifier.Trained{1}.predict( dataVectors.Z(k,:) );
        ws(2) = b(2);
        if strcmp( icPool.Z_v_C.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.Z_v_C.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Z_v_Y.classifier.Trained{1}.predict( dataVectors.Z(k,:) );
        ws(3) = b(2);
        if strcmp( icPool.Z_v_Y.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.Z_v_Y.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Z_v_O.classifier.Trained{1}.predict( dataVectors.Z(k,:) );
        ws(4) = b(2);
        if strcmp( icPool.Z_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.Z_v_O.classifier.NumTrainedPerFold;
        end
        
        %scores.Z(k) = min(ws);
        scores.Z(k,:) = ws;
    end
    
    scores.O = zeros(size( dataVectors.O, 1 ), 4);
    for k = 1 : size( dataVectors.O, 1 )
        ws = zeros(1, 4);

        [a, b] = icPool.B_v_O.classifier.Trained{1}.predict( dataVectors.O(k,:) );
        ws(1) = b(1);
        if strcmp( icPool.B_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(1) = ws(1) * icPool.B_v_O.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.C_v_O.classifier.Trained{1}.predict( dataVectors.O(k,:) );
        ws(2) = b(1);
        if strcmp( icPool.C_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(2) = ws(2) * icPool.C_v_O.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Y_v_O.classifier.Trained{1}.predict( dataVectors.O(k,:) );
        ws(3) = b(1);
        if strcmp( icPool.Y_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(3) = ws(3) * icPool.Y_v_O.classifier.NumTrainedPerFold;
        end
        [a, b] = icPool.Z_v_O.classifier.Trained{1}.predict( dataVectors.O(k,:) );
        ws(4) = b(1);
        if strcmp( icPool.Z_v_O.classifier.CrossValidatedModel, 'Bag' )
            ws(4) = ws(4) * icPool.Z_v_O.classifier.NumTrainedPerFold;
        end
        
        %scores.O(k) = min(ws);
        scores.O(k,:) = ws;
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
            
            if type == 1 || type == -1
                idx = find( ([dataPeaks.B.peakID] == peakID) & ([dataPeaks.B.spectrum] == spectra(s)) );
                if length(idx) == 1
                    weights(m) = scores.B(idx,1);
                elseif isempty(idx)
                    idx = find( ([dataPeaks.C.peakID] == peakID) & ([dataPeaks.C.spectrum] == spectra(s)) );
                    if length(idx) == 1
                        weights(m) = -scores.C(idx,1);
                    else
                        idx = find( ([dataPeaks.Y.peakID] == peakID) & ([dataPeaks.Y.spectrum] == spectra(s)) );
                        if length(idx) == 1
                            weights(m) = -scores.Y(idx,1);
                        else
                            idx = find( ([dataPeaks.Z.peakID] == peakID) & ([dataPeaks.Z.spectrum] == spectra(s)) );
                            if length(idx) == 1
                                weights(m) = -scores.Z(idx,1);
                            else
                                idx = find( ([dataPeaks.O.peakID] == peakID) & ([dataPeaks.O.spectrum] == spectra(s)) );
                                if length(idx) == 1
                                    weights(m) = -scores.O(idx,1);
                                else
                                    stop = 1;
                                end
                            end
                        end
                    end
                else
                    stop = 1;
                end
            elseif type == 2 || type == -2
                idx = find( ([dataPeaks.C.peakID] == peakID) & ([dataPeaks.C.spectrum] == spectra(s)) );
                if length(idx) == 1
                    weights(m) = scores.C(idx,1);
                elseif isempty(idx)
                    idx = find( ([dataPeaks.B.peakID] == peakID) & ([dataPeaks.B.spectrum] == spectra(s)) );
                    if length(idx) == 1
                        weights(m) = -scores.B(idx,2);
                    else
                        idx = find( ([dataPeaks.Y.peakID] == peakID) & ([dataPeaks.Y.spectrum] == spectra(s)) );
                        if length(idx) == 1
                            weights(m) = -scores.Y(idx,2);
                        else
                            idx = find( ([dataPeaks.Z.peakID] == peakID) & ([dataPeaks.Z.spectrum] == spectra(s)) );
                            if length(idx) == 1
                                weights(m) = -scores.Z(idx,2);
                            else
                                idx = find( ([dataPeaks.O.peakID] == peakID) & ([dataPeaks.O.spectrum] == spectra(s)) );
                                if length(idx) == 1
                                    weights(m) = -scores.O(idx,2);
                                else
                                    stop = 1;
                                end
                            end
                        end
                    end
                else
                    stop = 1;
                end
            end
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
        byPeaks, byIonClassifier, spectra(s).filename };
end
%disp( 'Finished' );
