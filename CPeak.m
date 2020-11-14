classdef CPeak < handle  % By Pengyu Hong @ Brandeis University
    properties
        mSpectrum = CSpectrum.empty(0,0); % For easy access.
        mMass = -1;  % Protonated mass
        mComplement = 0; % 0 - no complement, positive - existing complement, negative - computational complement
        mIsComplement = 0; % 0 - is not a complementary peak, 1 - otherwise
        mIntensity = 0.000001;  % normalized to sum up to 1
        mOriPeaks = CPeak.empty(0,0);
        
        mZscore = 0;
        type = 'O';
        RE = '';
        NRE = '';
        linkage = -1;
        mFeature = [];
        mFeatureIntensity = [];
        mNeighbour;
        mLinkageVector;
        ID;
        
        mInferredSuperSet = CTopologySuperSet.empty(0,0);
        mInferredFormulas = {};
        mInferredMasses = [];
        mInferredScores = [];
        
        mIonClassifierScores = []; % [B-ion score, C-ion score]
        
        mMassRange = [0, 0];
        mRawMZ = -1;
        mRawZ = -1;
        mComment = '';
        
        mType = CPeakType.Unknown;
    end
    
    methods
        function result = num(obj) % A CPeak handle may point to a CPeak array.
            result = length(obj);
        end
        
        function clear_inferred( obj )
            for k = 1 : length(obj)
                obj(k).mInferredSuperSet = CTopologySuperSet.empty(0,0);
                obj(k).mInferredFormulas = {};
                obj(k).mInferredMasses = [];
                obj(k).mInferredScores = [];
            end
        end
        
        function clear(obj)
            obj.mSpectrum = CSpectrum.empty(0,0);
            obj.mMass = -1;  % Protonated mass
            obj.mComplement = 0; % 0 - no complement, positive - existing complement, negative - computational complement
            obj.mIsComplement = 0; % 0 - is not a complementary peak, 1 - otherwise
            obj.mIntensity = 0.000001;  % normalized to sum up to 1
            obj.mOriPeaks = CPeak.empty(0,0);
            
            obj.mInferredSuperSet = CTopologySuperSet.empty(0,0);
            obj.mInferredFormulas = {};
            obj.mInferredMasses = [];
            obj.mInferredScores = [];
            
            obj.mIonClassifierScores = []; % [B-ion score, C-ion score]
            
            obj.mMassRange = [0, 0];
            obj.mRawMZ = -1;
            obj.mRawZ = -1;
            obj.mComment = '';
        end
        
        function copy( obj, aPeak )
            if isempty( aPeak )
                obj.clear();
            else
                fields = fieldnames( aPeak );
                for f = 1 : length( fields )
                    %if ( strcmp( fields{f}, 'mInferredSuperSet' ) )
                    %    if ~isempty( aPeak.(fields{f}) )
                    %        obj.(fields{f}) = aPeak.(fields{f}).copy;
                    %    else
                    %        obj.(fields{f}) = CTopologySuperSet.empty(0,0);
                    %    end
                    %else
                        obj.(fields{f}) = aPeak.(fields{f});
                    %end
                end
            end
        end
        
        function result = add_complementary_ions( obj, precursorMass, perMe, ionMass, minMassThreshold, massAccuracy )
            if isempty( obj )
                result = [];
                return;
            end
            if nargin < 2 || isempty( precursorMass )
                precursorMass = max( [obj.mMass] );
            end
            if nargin < 3 || isempty(perMe)
                perMe = 1; 
            end
            if nargin < 4 || isempty( ionMass )
                ionMass = CMass.Proton;
            end
            if nargin < 5 || isempty( minMassThreshold )
                minMassThreshold = 100;
            end
            if nargin < 6 || isempty( massAccuracy )
                massAccuracy = 5; % 5ppm
            end
            if perMe 
                minMassThreshold = max( minMassThreshold, CMonosaccharideSet.members(1).mPermethylated.mass - CMass.CH2 - CMass.H2O );
            else
                minMassThreshold = max( minMassThreshold, CMonosaccharideSet.members(1).mNative.mass - CMass.H2O );
            end
            
            for k = 1 : length(obj)
                if obj(k).mComplement > length(obj)
                    disp(k);
                end
            end
            
            masses = [obj.mMass];
            num = length( obj );
            complements = CPeak.empty(0, num);
            idx = 0;
            comD = ones(1, num) * 10000;
            for k = 1 : num
                temp = precursorMass - obj(k).mMass + ionMass;
                if temp < minMassThreshold
                    continue;
                end
                
                delta = (precursorMass + obj(k).mMass) * massAccuracy / 1000000; % 2019/9/11
                % delta = precursorMass * massAccuracy / 1000000; % Changed to the above on 2019/9/11
                [md, midx] = min( abs( masses - temp ) );
                if md > delta % not in the list
                    idx = idx + 1;
                    complements(idx) = CPeak();
                    complements(idx).mSpectrum = obj(k).mSpectrum;
                    complements(idx).mMass = temp;
                    complements(idx).mComplement = -k;
                    complements(idx).mIsComplement = 1;
                    complements(idx).mIntensity = obj(k).mIntensity;
                    complements(idx).mMassRange = temp + [-delta, delta];
                    if obj(k).mComplement == 0
                        obj(k).mComplement = -(length(obj) + idx);
                    end
                else
                    if comD(k) > md
                        %if obj(midx).mComplement > 0
                        %    obj(obj(midx).mComplement).mComplement = 0;
                        %end
                        obj(midx).mComplement = k;
                        obj(k).mComplement = midx;
                        comD(k) = md;
                        comD(midx) = md;
                    end
                    % complements(idx).mIsComplement = 0;
                end
            end
            
            if idx > 0
                result = [obj, complements(1:idx)];
                result = result.sort();
            else
                result = obj;
            end
            
            % A peak may have more than one complementary peaks
            % when the masses of the complementary peaks are
            % very close to each other. But only the closest
            % complementary peak is kept in record.
            % for k = 1 : length(result)
            %    if result(k).mComplement > 0
            %        if result( result(k).mComplement ).mComplement ~= k
                       % disp( [0, k, result(k).mComplement] ); 
                       % A peak may have more than one complementary peaks
                       % when the masses of the complementary peaks are
                       % very close to each other. But only the closest
                       % complementary peak is kept in record.
            %        end
            %    elseif result(k).mComplement < 0
            %        if result( -result(k).mComplement ).mComplement ~= -k
            %            disp( [1, k, result(k).mComplement] );
            %        end
            %    end
            % end
        end
        
        function [score, peakIds] = scoring( obj, spect, massAccuracy )
            if nargin < 3
                massAccuracy = 0.01; 
            end
            
            myMasses = unique( [obj.mMass] );
            if isa( spect, 'CSpectrum' )
                otherMasses = unique( [spect.mPeaks.mMass] );
            elseif isa( spect, 'CPeak' )
                otherMasses = unique( [spect.mMass] );
            else
                otherMasses = spect;
            end
            
            diff = ones(length(otherMasses), 1) * myMasses - otherMasses' * ones(1, length(myMasses));
            score = sum( abs( diff ) < massAccuracy ); % every column is a peak in obj.
            peakIds = find( score > 0 );
            score = length(peakIds);
        end
        
        function results = sort( obj, order ) % order = 'ascend'/'descend'
            if nargin < 2, order = 'ascend'; end
            [~, idx] = sort( [obj.mMass], order );
            results = obj(idx);
            for k = 1 : length(results)
                if results(k).mComplement > 0
                    results(k).mComplement = find( idx == results(k).mComplement );
                elseif results(k).mComplement < 0
                    results(k).mComplement = -find( idx == -results(k).mComplement );
                end
            end
        end
        
        function [peaks, idx] = find( obj, mz, err_threshold )
            if isempty( obj )
                peaks = []; idx = [];
            else
                % Updated on 2019/10/01
                if nargin < 3 || isempty(err_threshold) % Use mMassRange
                    idx = [];
                    for k = 1 : length( obj )
                        if (obj(k).mMassRange(1) <= mz) && (obj(k).mMassRange(2) >= mz)
                            idx(end+1) = k;
                        end
                    end
                    peaks = obj(idx);
                else
                    idx = find( abs([obj.mMass] - mz) < err_threshold );
                    peaks = obj(idx);
                end
                % Updated on 2019/10/01
            end            
        end
        
        function [peaks, idx] = find_rawmz( obj, mz, err_threshold )
            if isempty( obj )
                peaks = [];
                idx = [];
            else
                if nargin < 3 || isempty(err_threshold)
                    err_threshold = 0.01;
                end
                idx = find( abs([obj.mRawMZ] - mz) < err_threshold );
                peaks = obj(idx);
            end            
        end
        
        function result = merge_peaks( obj, threshold, massAccuracy )
            if nargin < 2, threshold = 0.001; end
            if nargin < 3, massAccuracy = 5; end % PPM

            old_mass = [obj.mMass]';
            old_intensity = [obj.mIntensity]';
            
            tree = linkage( old_mass, 'single' );
            c = cluster( tree, 'cutoff', threshold, 'criterion', 'distance' );
            newLen = max(c);
            new_mass = zeros(newLen, 1);
            new_intensity = zeros(newLen, 1);
            oldIdx = cell(newLen, 1);
            
            for k = 1 : newLen
                idx = find( c == k );
                tempM = old_mass( idx );
                tempI = old_intensity( idx );
                new_intensity(k) = sum( tempI );
                new_mass(k) = tempM' * tempI / new_intensity(k);
                oldIdx{k} = idx;
            end
            
            [temp, idx] = sortrows( [new_mass, new_intensity], 1 );
            oldIdx = oldIdx(idx);
            result = CPeak.create_peaks( temp(:,1), temp(:,2) );
            for k = 1 : length( result )
                result(k).mSpectrum = obj(1).mSpectrum;
                idx = oldIdx{k};
                if length( idx ) == 1
                    result(k).copy( obj(idx) ); % 2019/09/30 Replace the below code block.
                    
                    % Begin -- 2019/09/30 Replaced by the above 
                    % result(k).mRawMZ = obj(idx).mRawMZ;
                    % result(k).mRawZ = obj(idx).mRawZ;
                    % result(k).mComplement = obj(idx).mComplement; % this need to be redirected.
                    % result(k).mIsComplement = obj(idx).mIsComplement;
                    % result(k).mMassRange = obj(idx).mMassRange;
                    % End -- 2019/09/30 Replaced by the above 
                else
                    % Begin -- 2019/09/30 Replaced by the below 
                    % complementS = NaN;
                    % for m = 1 : length(idx)
                    %    if obj(idx(m)).mComplement > 0
                    %        complementS = idx(m);
                    %    elseif obj(idx(m)).mComplement < 0 && isnan(complementS)
                    %        complementS = idx(m);
                    %    elseif obj(idx(m)).mComplement == 0 && (isnan(complementS) || complementS < 0)
                    %        complementS = idx(m);
                    %    end
                    % end
                    % result(k).mComplement = obj(complementS).mComplement;
                    % result(k).mIsComplement = obj(complementS).mIsComplement;
                    % if result(k).mIsComplement && (result(k).mIsComplement > 0)
                    %    disp(['wrong ', num2str(k)]);
                    % end
                    % End -- 2019/09/30 Replaced by the below
                    
                    % Begin -- 2019/09/30 Replace the above
                    lowMass = inf;
                    highMass = -inf;
                    for m = 1 : length( idx )
                        oldP = obj(idx(m));
                        if oldP.mIsComplement == 0
                            if isempty( obj(idx(m)).mOriPeaks ) % is an original peak
                                result(k).mOriPeaks = [result(k).mOriPeaks, oldP];
                            else % is a merged peak.
                                result(k).mOriPeaks = [result(k).mOriPeaks, oldP.mOriPeaks];
                            end
                        end
                        if oldP.mMassRange(1) < lowMass
                            lowMass = oldP.mMassRange(1);
                        end
                        if oldP.mMassRange(2) > highMass
                            highMass = oldP.mMassRange(2);
                        end
                    end
                    if isempty( result(k).mOriPeaks )
                        result(k).mIsComplement = 1;
                        delta = massAccuracy / 1000000;
                        result(k).mMassRange = result(k).mMass * [1-delta, 1+delta];
                    else
                        result(k).mMassRange = [lowMass, highMass];
                    end
                    % End -- 2019/09/30 Replace the above
                end
            end
            
            % re-direct result(k).mComplement
            % Begin -- 2019/09/30 Replace the below
            masses = [result.mMass];
            precursorMass = result(k).mSpectrum.mPrecursor;
            delta = precursorMass * massAccuracy / 1000000;
            for k = 1 : length( result )
                temp = precursorMass - result(k).mMass + CMass.Proton;
                [minV, idx] = min( abs(masses - temp) );
                if idx > k
                    if minV <= delta
                        if result(k).mIsComplement == 1 || result(idx).mIsComplement == 1
                            result(k).mComplement = -idx;
                            result(idx).mComplement = -k;
                        else
                            result(k).mComplement = idx;
                            result(idx).mComplement = k;
                        end
                    elseif result(k).mIsComplement == 1
                        fprintf('Cannot find complementary after merging peaks: %d\n', k)
                    end
                end
            end
            % End -- 2019/09/30 Replace the below 

            % Begin -- 2019/09/30 Replaced by the above
            % masses = [result.mMass];
            % delta = result(k).mSpectrum.mPrecursor * massAccuracy / 1000000;
            % for k = 1 : length( result )
            %    if result(k).mComplement ~= 0
            %        temp = result(k).mSpectrum.mPrecursor - result(k).mMass + CMass.Proton;
            %        [minV, idx] = min( abs(masses - temp) );
            %        if minV <= delta
            %            result(k).mComplement = idx * sign(result(k).mComplement);
            %        else
            %            fprintf('Cannot find complementary after merging peaks: %d\n', k)
            %        end
            %    end
            % end
            % End -- 2019/09/30 Replaced by the above
        end
        
        function result = protonate( obj, considerOtherCharge, massAccuracy )
            if nargin < 2 || isempty( considerOtherCharge )
                considerOtherCharge = 0;
            end
            if nargin < 3 || isempty( massAccuracy )
                massAccuracy = 5; % PPM
            end
            
            % mz.1H = mz.Sodium * z - z * (CMass.Sodium - CMass.electron) + CMass.Proton
            newMasses = zeros( length(obj)*2, 2 ); 
            newMassIdx = 0;
            for k = 1 : length( obj )
                tmass = CPeak.protonate_raw_mz( obj(k).mRawMZ, obj(k).mRawZ, obj(k).mSpectrum.mMetal );
                obj(k).mMass = tmass(1);
                delta = obj(k).mMass * massAccuracy / 1000000;
                obj(k).mMassRange = obj(k).mMass + [-delta, delta];
                
                if considerOtherCharge && length(tmass) > 1
                    for m = 2 : length(tmass)
                        newMassIdx = newMassIdx + 1;
                        newMasses(newMassIdx, :) = [tmass(m), k];
                    end
                end
            end
            
            if newMassIdx > 0
                newPeaks = CPeak.create_peaks( newMasses(1:newMassIdx, 1) );
                for k = 1 : length( newPeaks )
                    newPeaks(k).mSpectrum = obj(1).mSpectrum;
                    newPeaks(k).mMass = newMasses(k, 1);
                    delta = newPeaks(k).mMass * massAccuracy / 1000000;
                    newPeaks(k).mMassRange = newPeaks(k).mMass + [-delta, delta];
                    sourcePeakIdx = newMasses(k, 2);
                    newPeaks(k).mIntensity = obj(sourcePeakIdx).mIntensity;
                    newPeaks(k).mRawMZ = obj(sourcePeakIdx).mRawMZ;
                    newPeaks(k).mRawZ = obj(sourcePeakIdx).mRawZ;
                end
                result = [obj, newPeaks];
                result = result.merge_peaks( 0.001, massAccuracy );
            else
                result = obj;
            end
            
            result = result.sort();
        end
    end
    
    methods (Static)
        function peaks = create_peaks( masses, intensities )
            num = length( masses );
            flag = 0;
            if nargin == 2 && length(intensities) == num
                flag = 1;
            end
            peaks = CPeak.empty( num, 0 );
            for k = 1 : num
                peaks(k) = CPeak;
                peaks(k).mMass = masses(k);
                if flag
                    peaks(k).mIntensity = intensities(k);
                end
            end
        end
        
        function mass = protonate_raw_mz( raw_mz, raw_z, metal )
            if raw_z == 0
                raw_z = 1;
            end
            temp = raw_mz * raw_z;
            metalMass = CMass.get_atom_mass( metal );
            mass = zeros(1, raw_z);
            
            % mass(1) = temp - raw_z * (metalMass - CMass.Electron) + CMass.Proton;
            for k = 1 : raw_z
                mass(k) = temp - k * (metalMass - CMass.Electron) - (raw_z - k) * CMass.Proton + CMass.Proton;
            end
            mass = unique( mass ); % if it is protonated, only one mass.
        end
    end
end