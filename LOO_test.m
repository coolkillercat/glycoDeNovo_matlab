function [CVResults, ICScores] = LOO_test( spectra, LOOResult )

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
                    w(1) = (2 * LOOResult.B_v_Y.LOO(1, idx) - 1) * LOOResult.B_v_Y.LOO(2, idx)';
                end
                idx = find( (LOOResult.B_v_Z.XInfo(:,1) == s) & (LOOResult.B_v_Z.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.B_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.B_v_Z.LOO(1, idx) - 1) * LOOResult.B_v_Z.LOO(2, idx)';
                end
                idx = find( (LOOResult.B_v_C.XInfo(:,1) == s) & (LOOResult.B_v_C.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.B_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.B_v_C.LOO(1, idx) - 1) * LOOResult.B_v_C.LOO(2, idx)';
                end
                idx = find( (LOOResult.B_v_O.XInfo(:,1) == s) & (LOOResult.B_v_O.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.B_v_O.LOO(2, idx);
                    w(4) = (2 * LOOResult.B_v_O.LOO(1, idx) - 1) * LOOResult.B_v_O.LOO(2, idx)';
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
                    w(1) = (2 * LOOResult.Y_v_O.LOO(1, idx) - 1) * LOOResult.Y_v_O.LOO(2, idx)';
                end
                idx = find( (LOOResult.Y_v_Z.XInfo(:,1) == s) & (LOOResult.Y_v_Z.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.Y_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.Y_v_Z.LOO(1, idx) - 1) * LOOResult.Y_v_Z.LOO(2, idx)';
                end
                idx = find( (LOOResult.Y_v_C.XInfo(:,1) == s) & (LOOResult.Y_v_C.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.Y_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.Y_v_C.LOO(1, idx) - 1) * LOOResult.Y_v_C.LOO(2, idx)';
                end
                idx = find( (LOOResult.Y_v_B.XInfo(:,1) == s) & (LOOResult.Y_v_B.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.Y_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.Y_v_B.LOO(1, idx) - 1) * LOOResult.Y_v_B.LOO(2, idx)';
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
                    w(1) = (2 * LOOResult.C_v_Y.LOO(1, idx) - 1) * LOOResult.C_v_Y.LOO(2, idx)';
                end
                idx = find( (LOOResult.C_v_Z.XInfo(:,1) == s) & (LOOResult.C_v_Z.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.C_v_Z.LOO(2, idx);
                    w(2) = (2 * LOOResult.C_v_Z.LOO(1, idx) - 1) * LOOResult.C_v_Z.LOO(2, idx)';
                end
                idx = find( (LOOResult.C_v_O.XInfo(:,1) == s) & (LOOResult.C_v_O.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.C_v_O.LOO(2, idx);
                    w(3) = (2 * LOOResult.C_v_O.LOO(1, idx) - 1) * LOOResult.C_v_O.LOO(2, idx)';
                end
                idx = find( (LOOResult.C_v_B.XInfo(:,1) == s) & (LOOResult.C_v_B.XInfo(:,2) == peakID) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.C_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.C_v_B.LOO(1, idx) - 1) * LOOResult.C_v_B.LOO(2, idx)';
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
                    w(1) = (2 * LOOResult.Z_v_O.LOO(1, idx) - 1) * LOOResult.Z_v_O.LOO(2, idx)';
                end
                idx = find( (LOOResult.Z_v_Y.XInfo(:,1) == s) & (LOOResult.Z_v_Y.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(2) = 0;
                else
                    % w(2) = LOOResult.Z_v_Y.LOO(2, idx);
                    w(2) = (2 * LOOResult.Z_v_Y.LOO(1, idx) - 1) * LOOResult.Z_v_Y.LOO(2, idx)';
                end
                idx = find( (LOOResult.Z_v_C.XInfo(:,1) == s) & (LOOResult.Z_v_C.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(3) = 0;
                else
                    % w(3) = LOOResult.Z_v_C.LOO(2, idx);
                    w(3) = (2 * LOOResult.Z_v_C.LOO(1, idx) - 1) * LOOResult.Z_v_C.LOO(2, idx)';
                end
                idx = find( (LOOResult.Z_v_B.XInfo(:,1) == s) & (LOOResult.Z_v_B.XInfo(:,2) == -spectra(s).mPeaks(peakID).mComplement) );
                if isempty(idx)
                    missed = missed + 1;
                    w(4) = 0;
                else
                    % w(4) = LOOResult.Z_v_B.LOO(2, idx);
                    w(4) = (2 * LOOResult.Z_v_B.LOO(1, idx) - 1) * LOOResult.Z_v_B.LOO(2, idx)';
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