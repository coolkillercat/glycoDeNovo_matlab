classdef CSpectrum < handle % By Pengyu Hong @ Brandeis University
    properties
        mPrecursor = -1 % mz of the protonated precursor
        mPeaks = CPeak.empty(0,0);
        
        mExperimentMethod = '';
        mMetal = '';
        mNLinked = 0;
        mPermethylated = 0;
        mProtonated = 0;
        mMassAccuracy = 5; %PPM

        mReducingEndModification = '';
        mComposition = zeros(1, CMonosaccharideSet.cNumberMonosaccharideClasses );
        comment = '';
        filename = '';
        
        mTopologies = CTopology.empty(0);
    end
    
    methods
        function result = num( obj )
            result = length(obj.mPeaks);
        end
        
        function clear_inferred( obj )
            for k = 1 : length(obj.mPeaks)
                obj.mPeaks(k).clear_inferred();
            end
        end
        
        function copy( obj, spec )
            if ~isempty( spec )
                fields = fieldnames( spec(1) );
                for f = 1 : length( fields )
                    if strcmp( fields{f}, 'mPeaks' )
                        num = length(spec.mPeaks);
                        obj.mPeaks = CPeak.empty(0,num);
                        for k = 1 : num
                            obj.mPeaks(k) = CPeak;
                            obj.mPeaks(k).copy( spec.mPeaks(k) );
                        end
                    else
                        obj.(fields{f}) = spec.(fields{f});
                    end
                end
            else
                obj.mPrecursor = -1;
                obj.mPeaks = CPeak.empty(0,0);
                
                obj.mExperimentMethod = '';
                obj.mMetal = '';
                obj.mNLinked = 0;
                obj.mPermethylated = 0;
                obj.mProtonated = 0;
                obj.mMassAccuracy = 5; %PPM
                
                obj.mReducingEndModification = '';
                obj.comment = '';
                obj.filename = '';
            end
        end
        
        function add_complementary_ions( obj, ionMass, minMassThreshold )
            if nargin < 2 || isempty( ionMass )
                ionMass = CMass.Proton;
            end
            if nargin < 3 || isempty( minMassThreshold )
                minMassThreshold = 100;
            end
            
            obj.mPeaks = obj.mPeaks.add_complementary_ions( obj.mPrecursor, obj.mPermethylated, ionMass, minMassThreshold );
        end
        
        function masses = get_masses( obj )
            if isempty( obj.mPeaks )
                masses = [];
            else
                masses = [obj.mPeaks.mMass];
            end
        end
        
        function [score, peakIds] = scoring( obj, spect, massAccuracy )
            if nargin < 3
                massAccuracy = 0.01; 
            end
            myMasses = unique( [obj.mPeaks.mMass] );
            if isa( spect, 'CSpectrum' )
                otherMasses = unique( [spect.mPeaks.mMass] );
            elseif isa( spect, 'CPeak' )
                otherMasses = unique( [spect.mPeaks.mMass] );
            else
                otherMasses = spect;
            end
            diff = ones(length(otherMasses), 1) * myMasses - otherMasses' * ones(1, length(myMasses));
            score = sum( abs( diff ) < massAccuracy ); % every column is a peak in obj.
            peakIds = find( score > 0 );
            score = length(peakIds);
        end
        
        function sort_peaks( obj, order ) % order = 'ascend'/'descend'
            if nargin < 2, order = 'ascend'; end
            obj.mPeaks = obj.mPeaks.sort( order );
        end
        
        function [peaks, idx] = find_peaks( obj, mass, massAccuracy )
            if nargin < 3 || isempty( massAccuracy )
                massAccuracy = 0.01;
            end
            if isempty( obj.mPeaks )
                peaks = [];
                idx = [];
            else
                [peaks, idx] = obj.mPeaks.find( mass, massAccuracy );
            end
        end
        
        function [peaks, idx] = find_peaks_rawmz( obj, mz, massThreshold )
            if nargin < 3 || isempty( massThreshold )
                massThreshold = 0.005;
            end
            if isempty( obj.mPeaks )
                peaks = [];
                idx = [];
            else
                num = length(obj.mPeaks);
                flag = zeros(1, num);
                for k = 1 : num
                    aPeak = obj.mPeaks(k);
                    if ~isempty( aPeak.mOriPeaks )
                        for oPeak = aPeak.mOriPeaks
                            if abs( oPeak.mRawMZ - mz ) < massThreshold
                                flag( k ) = 1;
                                break;
                            end
                        end
                    elseif abs( obj.mPeaks(k).mRawMZ - mz ) < massThreshold
                        flag( k ) = 1;
                    end
                end
                idx = find( flag );
                peaks = obj.mPeaks(idx);
                % [peaks, idx] = obj.mPeaks.find_rawmz( mz, massThreshold );
            end
        end
        
        function merge_peaks( obj, threshold )
            obj.mPeaks = obj.mPeaks.merge_peaks( threshold );
        end
        
        function protonate( obj, considerOtherCharge )
            if obj.mProtonated
                return;
            end
            if nargin < 2
                considerOtherCharge = 1; 
            end
            obj.mPeaks = obj.mPeaks.protonate( considerOtherCharge, obj.mMassAccuracy );
            obj.mPeaks = obj.mPeaks( [obj.mPeaks.mMass] <= obj.mPrecursor * (1 + obj.mMassAccuracy/1000000) );
            obj.mProtonated = 1;
        end
        
        function result = replicate( obj )
            result = CSpectrum;
            fields = fieldnames( obj );
            for f = 1 : length( fields )
                if strcmp( fields{f}, 'mPeaks' )
                    num = length( obj.mPeaks );
                    result.mPeaks = CPeak.empty(0, num);
                    for k = 1 : num
                        result.mPeaks(k) = CPeak;
                        result.mPeaks(k).copy( obj.mPeaks(k) );
                    end
                else
                    result.(fields{f}) = obj.(fields{f});
                end
            end
        end
        
        function rankInfo = get_rank(this, glycanName)
            rankInfo.rankBySupportingPeak = 0;
            rankInfo.corankBySupportingPeak = 0;
            rankInfo.rankByIonclassifier = 0;
            rankInfo.corankByIonclassifier = 0;
            rankInfo.idx = 0;
            candidates = {};
            scores = [];
            for TSS = this.mPeaks(end).mInferredSuperSet
                for TP = TSS.mTopologies
                    candidates{end+1} = TP.mFormula;
                    scores(end+1,:) = [length(TP.mSupportPeaks), TP.mScore];
                end
            end
            idx = find(strcmp( candidates, glycanName ));
            if isempty(idx)
                candidates = {};
                for TSS = this.mPeaks(end).mInferredSuperSet
                    for TP = TSS.mTopologies
                        g = CGlycan;
                        g.parse( TP.mFormula );
                        g.order_branch_by_name();
                        candidates{end+1} = g.mInferredFormula;
                    end
                end
                idx = find(strcmp( candidates, glycanName ));
            end
            if ~isempty(idx)
                rankInfo.rankBySupportingPeak = sum( scores(:,1) > scores(idx,1) ) + 1;
                rankInfo.corankBySupportingPeak = sum( scores(:,1) == scores(idx,1) );
                rankInfo.rankByIonclassifier = sum( scores(:,2) > scores(idx,2) ) + 1;
                rankInfo.corankByIonclassifier = sum( scores(:,2) == scores(idx,2) );
                rankInfo.idx = idx;
                temp = unique(scores(:,2));
                if length(temp) > 1
                    rankInfo.topIonCalssifierScores = temp(end-1:end);
                else
                    rankInfo.topIonCalssifierScores = [NaN, temp];
                end
            end
        end
        
        function save_reconstruction( obj, filename, ion_filter, checkMinus2H, allowGap )
            if nargin < 2 || isempty( filename )
                disp( 'Please specify a filename.' ); return
            end
            if nargin < 3 || isempty( ion_filter )
                ion_filter = {'Y', 'B', 'X', 'A', 'Z', 'C', 'T'};
            end
            if nargin < 4, checkMinus2H = []; end
            if nargin < 5, allowGap = []; end
            
            fid = fopen( filename, 'w' );
            flag = 0;
            if ~isempty( obj.comment )
                fprintf( fid, '# %s\n', obj.comment );
                flag = 1;
            end
            if ~isempty( checkMinus2H )
                fprintf( fid, '# Check -2H = %d\n', checkMinus2H );
                flag = 1;
            end
            if ~isempty( allowGap )
                fprintf( fid, '# Allow gap = %d\n', allowGap );
                flag = 1;
            end
            if any( obj.mComposition )
                fprintf( fid, '# Composition:' );
                for k = 1 : length( obj.mComposition )
                    if obj.mComposition(k) > 0
                        fprintf( fid, ' [%s: %d]', CMonosaccharideSet.cMonoClasses{k}, obj.mComposition(k) );
                    end
                end
                fprintf( fid, '\n' );
                flag = 1;
            end
            if flag, fprintf( fid, '\n' ); end
            
            for k = 1 : length( obj.mPeaks )
                aPeak = obj.mPeaks(k);
                flag = 0;
                for m = 1 : length( aPeak.mInferredSuperSet )
                    flag = sum( strcmp( ion_filter, aPeak.mInferredSuperSet(m).mType ) ) && ...
                        ~isempty( aPeak.mInferredFormulas );
                    if flag, break; end
                end
                if flag == 0, continue; end
                
                if aPeak.mIsComplement
                    if aPeak.mComplement ~= 0
                        cIdx = -aPeak.mComplement;
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mComplement, aPeak.mMass, obj.mPeaks(cIdx).mRawMZ, obj.mPeaks(cIdx).mRawZ, ...
                                aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mComplement, aPeak.mMass, obj.mPeaks(cIdx).mRawMZ, obj.mPeaks(cIdx).mRawZ, ...
                                aPeak.mIntensity, aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                    else
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d [IsCom] mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d [IsCom] mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, aPeak.mIntensity, ...
                                aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                    end
                elseif aPeak.mComplement
                    if ~isempty( obj.mPeaks(k).mOriPeaks )
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mComplement, aPeak.mMass, aPeak.mOriPeaks(1).mRawMZ, ...
                                aPeak.mOriPeaks(1).mRawZ, aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mComplement, aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, ...
                                aPeak.mIntensity, aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                        
                    else
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mComplement, aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mComplement, aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, ...
                                aPeak.mIntensity, aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                    end
                else
                    if ~isempty( aPeak.mOriPeaks )
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mMass, aPeak.mOriPeaks(1).mRawMZ, aPeak.mOriPeaks(1).mRawZ, aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mMass, aPeak.mOriPeaks(1).mRawMZ, aPeak.mOriPeaks(1).mRawZ, aPeak.mIntensity, ...
                                aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                    else
                        if isempty( aPeak.mIonClassifierScores )
                            fprintf( fid, '@ Peak %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                                aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, aPeak.mIntensity);
                        else
                            fprintf( fid, '@ Peak %d mass %f rawMZ %f rawZ %d intensity %d. IC = [%.2f, %.2f]\n', k, ...
                                aPeak.mMass, aPeak.mRawMZ, aPeak.mRawZ, aPeak.mIntensity, ...
                                aPeak.mIonClassifierScores(1), aPeak.mIonClassifierScores(2));
                        end
                    end
                end
                
                if k == length(obj.mPeaks) && ~isempty( obj.mTopologies )
                    for aTP = obj.mTopologies
                        fprintf( fid, '** T: %s [Peaks (%d, %.2f):', ...
                            aTP.mFormula, length( aTP.mSupportPeaks ), aTP.mScore );
                        fprintf( fid, ' %d', aTP.mSupportPeaks );
                        fprintf( fid, ']\n' );
                    end
                else
                    for p = 1 : length( obj.mPeaks(k).mInferredSuperSet )
                        pidx =  obj.mPeaks(k).mInferredSuperSet(p).mTargetPeaks(1,:) == k ;
                        tpType = obj.mPeaks(k).mInferredSuperSet(p).mTargetPeaks(2, pidx);
                        minusH = 0;
                        switch tpType
                            case CPeakType.B
                                tpType = 'B';
                            case CPeakType.B_H
                                tpType = 'B';
                                minusH = 1;
                            case CPeakType.B_2H
                                tpType = 'B';
                                minusH = 2;
                            case CPeakType.C
                                tpType = 'C';
                            case CPeakType.C_H
                                tpType = 'C';
                                minusH = 1;
                            case CPeakType.C_2H
                                tpType = 'C';
                                minusH = 2;
                            case CPeakType.T
                                tpType = 'T';
                            case CPeakType.T_H
                                tpType = 'T';
                                minusH = 1;
                            case CPeakType.T_2H
                                tpType = 'T';
                                minusH = 2;
                        end
                        
                        for aTP = obj.mPeaks(k).mInferredSuperSet(p).mTopologies
                            if minusH
                                fprintf( fid, '** %s-%dH: %s [Peaks (%d, %.2f):', ...
                                    tpType, minusH, aTP.mFormula, ...
                                    length( aTP.mSupportPeaks ), aTP.mScore );
                            else
                                fprintf( fid, '** %s: %s [Peaks (%d, %.2f):', ...
                                    tpType, aTP.mFormula, ...
                                    length( aTP.mSupportPeaks ), aTP.mScore );
                            end
                            fprintf( fid, ' %d', aTP.mSupportPeaks );
                            fprintf( fid, ']\n' );
                        end
                    end
                end
                fprintf( fid, '\n' );
            end
            fclose( fid );
        end
        
        function save_reconstruction_full( obj, filename, checkMinus2H, allowGap )
            if nargin < 2 || isempty( filename )
                disp( 'Please specify a filename.' ); return
            end
            if nargin < 3, checkMinus2H = []; end
            if nargin < 4, allowGap = []; end
            
            fid = fopen( filename, 'w' );
            flag = 0;
            if ~isempty( obj.comment )
                fprintf( fid, '# %s\n', obj.comment );
                flag = 1;
            end
            if ~isempty( checkMinus2H )
                fprintf( fid, '# Check -2H = %d\n', checkMinus2H );
                flag = 1;
            end
            if ~isempty( allowGap )
                fprintf( fid, '# Allow gap = %d\n', allowGap );
                flag = 1;
            end
            if any( obj.mComposition )
                fprintf( fid, '# Composition:' );
                for k = 1 : length( obj.mComposition )
                    if obj.mComposition(k) > 0
                        fprintf( fid, ' [%s - %d]', CMass2Composition.mMonoClasses{k}, obj.mComposition(k) );
                    end
                end
            end
            if flag, fprintf( fid, '\n' ); end
            
            for k = 1 : length( obj.mPeaks )
                if obj.mPeaks(k).mIsComplement
                    cIdx = -obj.mPeaks(k).mComplement;
                    %{
                    fprintf( fid, '@ Peak %d (~ %d): mass %f [%f(%d)], intensity %d, IC = [%.2f, %.2f]\n', k, ...
                        obj.mPeaks(k).mComplement, obj.mPeaks(k).mMass, obj.mPeaks(cIdx).mRawMZ, obj.mPeaks(cIdx).mRawZ, ...
                        obj.mPeaks(k).mIntensity, obj.mPeaks(k).mIonClassifierScores(1), obj.mPeaks(k).mIonClassifierScores(2) );
                    %}
                    fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                        obj.mPeaks(k).mComplement, obj.mPeaks(k).mMass, obj.mPeaks(cIdx).mRawMZ, obj.mPeaks(cIdx).mRawZ, ...
                        obj.mPeaks(k).mIntensity);
                elseif obj.mPeaks(k).mComplement
                    %{
                    fprintf( fid, '@ Peak %d (~ %d): mass %f [%f(%d)], intensity %d, IC = [%.2f, %.2f]\n', k, ...
                        obj.mPeaks(k).mComplement, obj.mPeaks(k).mMass, obj.mPeaks(k).mRawMZ, obj.mPeaks(k).mRawZ, obj.mPeaks(k).mIntensity, ...
                        obj.mPeaks(k).mIonClassifierScores(1), obj.mPeaks(k).mIonClassifierScores(2) );
                    %}
                    fprintf( fid, '@ Peak %d complement %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                        obj.mPeaks(k).mComplement, obj.mPeaks(k).mMass, obj.mPeaks(k).mRawMZ, obj.mPeaks(k).mRawZ, obj.mPeaks(k).mIntensity);
                else
                    %{
                    fprintf( fid, '@ Peak %d : mass %f [%f(%d)], intensity %d, IC = [%.2f, %.2f]\n', k, ...
                        obj.mPeaks(k).mMass, obj.mPeaks(k).mRawMZ, obj.mPeaks(k).mRawZ, obj.mPeaks(k).mIntensity, ...
                        obj.mPeaks(k).mIonClassifierScores(1), obj.mPeaks(k).mIonClassifierScores(2) );
                    %}
                    fprintf( fid, '@ Peak %d mass %f rawMZ %f rawZ %d intensity %d\n', k, ...
                        obj.mPeaks(k).mMass, obj.mPeaks(k).mRawMZ, obj.mPeaks(k).mRawZ, obj.mPeaks(k).mIntensity);
                end
                
                fprintf( fid, '\n' );
            end
            fprintf( fid, '\n\n' );
            for k = 1 : length( obj.mPeaks )
                if isempty( obj.mPeaks(k).mInferredSuperSet ) || isempty(obj.mPeaks(k).mInferredSuperSet(1).mTopologies)
                    continue;
                end
                for p = 1 : length( obj.mPeaks(k).mInferredSuperSet )
                    pidx =  obj.mPeaks(k).mInferredSuperSet(p).mTargetPeaks(1,:) == k ;
                    switch obj.mPeaks(k).mInferredSuperSet(p).mTargetPeaks(2, pidx);
                        case CPeakType.B
                            atype = 'B';
                        case CPeakType.B_H
                            atype = 'B-1H';
                        case CPeakType.B_2H
                            atype = 'B-2H';
                        case CPeakType.C
                            atype = 'C';
                        case CPeakType.C_H
                            atype = 'C-1H';
                        case CPeakType.C_2H
                            atype = 'C-2H';
                        otherwise
                            atype = 'T';
                    end
                    for aTP = obj.mPeaks(k).mInferredSuperSet(p).mTopologies
                        fprintf( fid, '** Peak %d type %s Formula %s \\Formula Peaks ', k, atype, aTP.mFormula );
                        fprintf( fid, ' %d', aTP.mSupportPeaks );
                        fprintf( fid, ' \\Peaks Scores %d %.2f\n', length( aTP.mSupportPeaks ), aTP.mScore );
                    end
                end
            end
            fclose( fid );
            disp('reconstruct_full complete');
        end
        
        function sort_topologies_by_score(this)
            if isempty( this.mPeaks(end).mInferredFormulas )
                return
            end
            
            this.mTopologies = CTopology.empty(0);
            scores = [];
            for aTTS = this.mPeaks(end).mInferredSuperSet
                for aTP = aTTS.mTopologies
                    this.mTopologies(end+1) = aTP;
                    scores(end+1) = aTP.mScore;
                end
            end
            [~, idx] = sort( scores, 'descend' );
            this.mTopologies = this.mTopologies(idx);
        end
        
        function missing = compare_to_GWB( obj, filename ) % this function needs updates
            fid = fopen( filename, 'r' );
            data = textscan( fid, '%s%s%f%s%f%f' );
            fclose(fid);
            data = data{3};
            mymass = [obj.mass];
            
            flag = ones(1, length(data));
            for k = 1 : length( data )
                flag(k) = min( abs( data(k) - mymass ) ) > 0.0001;
            end
            missing.not_in_obj = data(flag >0);
            
            flag = ones(1, length(mymass));
            for k = 1 : length( mymass )
                flag(k) = min( abs( mymass(k) - data ) ) > 0.0001;
            end
            missing.not_in_gwb = mymass(flag >0);
        end
        
        function visualize(obj)
            masses = [obj.mPeaks.mMass];
            intensities = [obj.mPeaks.mIntensity];
            intensities(end) = max(intensities(1:end-1)) + 100;
            minMass = obj.mPeaks(1).mMass - 10;
            maxMass = obj.mPeaks(end).mMass + 10;
            intensities = log(intensities);
            bar( masses, intensities );
            axis( [minMass, maxMass, 0, intensities(end)] );
        end
        
        function G = visualize_reconstruction_graph(obj)
            G = digraph();
            peakSet = [];
            for TSS = obj.mPeaks(end).mInferredSuperSet
                for TP = TSS.mTopologies
                    peakSet = [peakSet, TP.mSupportPeaks];
                end
            end
            peakSet = unique( peakSet );
            
            vidx = 0;
            for k = peakSet
                % target = ['p', num2str(k)];
                target = num2str(k);
                for TSS = obj.mPeaks(k).mInferredSuperSet
                    for TPS = TSS.mTopologySets
                        vidx = vidx + 1;
                        vnode = ['v', num2str(vidx)];
                        % G = G.addedge( vnode, target );
                        for sIdx = 1 : 4
                            if isempty( TPS.mSources{sIdx} )
                                break;
                            end
                            tidx = 1;
                            while tidx < numel(TPS.mSources{sIdx}.mTargetPeaks)
                                % source = ['p', num2str(TPS.mSources{sIdx}.mTargetPeaks(tidx))];
                                source = num2str(TPS.mSources{sIdx}.mTargetPeaks(tidx));
                                tidx = tidx + 2;
                                if isempty( G.Nodes ) || ~G.findnode(target) || ~G.findnode(source) || ~G.findedge( {target}, {source} )
                                    G = G.addedge( target, source );
                                end
%                                 if isempty( G.Nodes ) || ~G.findnode(target) || ~G.findnode(vnode) || ~G.findedge( {target}, {vnode} )
%                                     G = G.addedge( target, vnode );
%                                 end
                            end
                        end
                    end
                end
            end
            p = plot(G);
            p.MarkerSize = 7;
        end
    end
    
    methods (Static)
        function result = create_spectrum( masses, intensities )
            result = CSpectrum;
            result.mPeaks = CPeak.create_peaks( masses, intensities );
        end
        
        function result = load( filename )
            chargeMetal = '';
            expMethod = '';
            per = 0;
            nlinked = 0;
            precursorMZ = -1;
            massAccuracy = 5; %PPM
            
            fid = fopen( filename, 'r' );
            if fid == -1, disp( ['no such file: ', filename] ); return; end
            
            result = CSpectrum;
            result.filename = filename;

            aline = fgetl( fid );
            firstLine = 1;
            while str_startswith( aline, '#' ) || ~isempty(strtrim(aline))
                if str_startswith( aline, '# Metal:' )
                    chargeMetal = strtrim( aline(9:end) );    
                    if strcmp( chargeMetal, 'H' )
                        chargeMetal = 'Proton';
                    end
                elseif str_startswith( aline, '# Method:' )
                    expMethod = strtrim( aline(10:end) );
                elseif str_startswith( aline, '# PPM:' )
                    massAccuracy = str2double( aline(7:end) );
                elseif str_startswith( aline, '# O18' )
                    result.mReducingEndModification = 'O18';
                elseif str_startswith( aline, '# Permethylated' )
                    per = 1;
                elseif strcmp( aline, '# PA' ) || strcmpi( aline, '# Aminopyridine' )
                    result.mReducingEndModification = 'Aminopyridine';
                elseif strcmp( aline, '# PRAGS' )
                    result.mReducingEndModification = 'PRAGS';
                elseif str_startswith( aline, '# Reduced' )
                    result.mReducingEndModification = 'Reduced';
                elseif str_startswith( aline, '# REM_' )
                    result.mReducingEndModification = aline(3:end);
                elseif str_startswith( aline, '# NLinked' )
                    nlinked = 1;
                elseif str_startswith( aline, '# Deuterium' ) || str_startswith( aline, '# D-Red' )
                    result.mReducingEndModification = 'Deuterium';
                elseif str_startswith( aline, '# Precursor:' )
                    temp = strtrim( aline(13:end) );
                    fields = strsplit( temp, ';' );
                    if strcmp( chargeMetal, 'Proton' )
                        temp = 'H+';
                    else
                        temp = [chargeMetal, '+'];
                    end
                    for k = 1 : length(fields)
                        fields{k} = strtrim( fields{k} );
                        idx = strfind( fields{k}, temp );
                        if ~isempty( idx )
                            precursorMZ = str2double( fields{k}(length(temp)+idx:end) );
                            if idx == 1
                                num = 1;
                            else
                                num = str2double( fields{k}(1:idx-1) );
                            end
                            precursorMZ = precursorMZ * num - num * (CMass.get_atom_mass( chargeMetal ) - CMass.Electron) + CMass.Proton;
                            break;
                        end
                    end
                elseif firstLine
                    result.comment = strtrim( aline(2:end) );
                    firstLine = 0;
                end
                
                aline = fgetl( fid );
                if ~ischar(aline), break, end
            end
            while ~str_startswith( aline, 'm/z' )
                aline = fgetl( fid );
                if ~ischar(aline), break, end
            end
            
            data = textscan( fid, '%f%s%f%*s' );
            fclose( fid );
            
            temp = zeros(length(data{2}), 1);
            for k = 1 : length(data{2})
                temp(k) = str2double( data{2}{k}(1) );
            end
            data = [data{1}, temp, data{3}];
            data = sortrows( data, 1 );
            
            result.mMetal = chargeMetal;
            result.mPermethylated = per;
            result.mExperimentMethod = expMethod;
            result.mNLinked = nlinked;
            result.mMassAccuracy = massAccuracy;
            
            result.mPeaks = CPeak.empty(size(data, 1), 0);
            for k = 1 : size(data, 1)
                result.mPeaks(k) = CPeak();
                result.mPeaks(k).mSpectrum = result;
                result.mPeaks(k).mRawMZ = data(k,1);
                result.mPeaks(k).mRawZ = data(k,2);
                result.mPeaks(k).mIntensity = data(k,3);
                if result.mProtonated
                    result.mPeaks(k).mMass = data(k,1);
                    delta = result.mPeaks(k).mMass * massAccuracy / 1000000;
                    result.mPeaks(k).mMassRange = result.mPeaks(k).mMass + [-delta, delta];
                end
            end
            
            if precursorMZ > 0
                result.mPrecursor = precursorMZ;
            elseif ~isempty( result.mPeaks )
                result.mPrecursor = max( [result.mPeaks.mMass] );
            end
        end
        
        function result = load_ann( filename )
            chargeMetal = '';
            expMethod = '';
            per = 0;
            nlinked = 0;
            precursorMZ = -1;
            massAccuracy = 5; %PPM
            
            fid = fopen( filename, 'r' );
            if fid == -1, disp( ['no such file: ', filename] ); return; end
            
            result = CSpectrum;
            result.filename = filename;

            aline = fgetl( fid );
            firstLine = 1;
            while str_startswith( aline, '#' ) || ~isempty(strtrim(aline))
                if str_startswith( aline, '# Metal:' )
                    chargeMetal = strtrim( aline(9:end) );    
                    if strcmp( chargeMetal, 'H' )
                        chargeMetal = 'Proton';
                    end
                elseif str_startswith( aline, '# Method:' )
                    expMethod = strtrim( aline(10:end) );
                elseif str_startswith( aline, '# PPM:' )
                    massAccuracy = str2double( aline(7:end) );
                elseif str_startswith( aline, '# O18' )
                    result.mReducingEndModification = 'O18';
                elseif str_startswith( aline, '# Permethylated' )
                    per = 1;
                elseif strcmp( aline, '# PA' ) || strcmpi( aline, '# Aminopyridine' )
                    result.mReducingEndModification = 'Aminopyridine';
                elseif strcmp( aline, '# PRAGS' )
                    result.mReducingEndModification = 'PRAGS';
                elseif str_startswith( aline, '# Reduced' )
                    result.mReducingEndModification = 'Reduced';
                elseif str_startswith( aline, '# REM_' )
                    result.mReducingEndModification = aline(3:end);
                elseif str_startswith( aline, '# NLinked' )
                    nlinked = 1;
                elseif str_startswith( aline, '# Deuterium' ) || str_startswith( aline, '# D-Red' )
                    result.mReducingEndModification = 'Deuterium';
                elseif str_startswith( aline, '# Precursor:' )
                    temp = strtrim( aline(13:end) );
                    fields = strsplit( temp, ';' );
                    if strcmp( chargeMetal, 'Proton' )
                        temp = 'H+';
                    else
                        temp = [chargeMetal, '+'];
                    end
                    for k = 1 : length(fields)
                        fields{k} = strtrim( fields{k} );
                        idx = strfind( fields{k}, temp );
                        if ~isempty( idx )
                            precursorMZ = str2double( fields{k}(length(temp)+idx:end) );
                            if idx == 1
                                num = 1;
                            else
                                num = str2double( fields{k}(1:idx-1) );
                            end
                            precursorMZ = precursorMZ * num - num * (CMass.get_atom_mass( chargeMetal ) - CMass.Electron) + CMass.Proton;
                            break;
                        end
                    end
                elseif firstLine
                    result.comment = strtrim( aline(2:end) );
                    firstLine = 0;
                end
                
                aline = fgetl( fid );
                if ~ischar(aline), break, end
            end
            while ~str_startswith( aline, 'm/z' )
                aline = fgetl( fid );
                if ~ischar(aline), break, end
            end
            
            result.mMetal = chargeMetal;
            result.mPermethylated = per;
            result.mExperimentMethod = expMethod;
            result.mNLinked = nlinked;
            result.mMassAccuracy = massAccuracy;
            
            data = fgetl(fid);
            x = 0;
            while ischar(data)
                x = x+1;
                datas = strsplit(data, '\t');
                data = fgetl(fid);
                currPeak = CPeak();
                currPeak.mSpectrum = result;
                currPeak.mRawMZ = str2double(datas{1});
                currPeak.mMass = currPeak.mRawMZ;
                currPeak.mRawZ = str2double(datas{2});
                currPeak.mZscore = str2double(datas{3});
                if result.mProtonated
                    result.mPeaks(k).mMass = data(k,1);
                    delta = result.mPeaks(k).mMass * massAccuracy / 1000000;
                    result.mPeaks(k).mMassRange = result.mPeaks(k).mMass + [-delta, delta];
                end
                if (strcmp(datas{4}, 'com') == 1)
                    currPeak.mIsComplement = 1;
                end
                if (length(datas) > 5)
                    if (strcmp(datas{4}, 'com') == 1)
                        start = 5;
                    else
                        start = 4;
                    end
                    currPeak.type = datas{start};
                    currPeak.RE = datas{start+1};
                    currPeak.NRE = datas{start+2};
                    if start + 3 == length(datas)
                        currPeak.linkage = 0;
                        currPeak.mComment = datas{start+3};
                    else
                        currPeak.linkage = str2double(datas{start+3});
                        currPeak.mComment = datas{start+4};
                    end
                end
                result.mPeaks(end+1) = currPeak;
                currPeak.ID = length(result.mPeaks);
            end
            
            if precursorMZ > 0
                result.mPrecursor = precursorMZ;
            elseif ~isempty( result.mPeaks )
                result.mPrecursor = max( [result.mPeaks.mMass] );
            end
            fclose(fid);
        end
    end
end