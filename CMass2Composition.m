classdef CMass2Composition < handle % By Pengyu Hong @ Brandeis University
% Use this after protonation but before anything else.
    properties (Constant)
        mMonoMaxNums = [4, 4, 9, 5, 8, 5, 7, 7];
    end
    
    properties (Access = protected)
        mPermethylated;
        mReducingEndModification = [];
        mREMMassCompensation = 0;
    end
    
    properties
        mCheckMinus2H;
        mMassAccuracyPPM;
        mMonoMasses;
        mPrecursorMassTolerance;
        mMasses;
        mCompositions;
    end
    
    methods
        function this = CMass2Composition()
            this.mCheckMinus2H = 0;
            this.mPrecursorMassTolerance = 0.01;
            this.mMassAccuracyPPM = 5;
            this.set_permethylation( 1 );
            this.mReducingEndModification = '';
            this.mREMMassCompensation = 0;
            this.mMasses = [];
            this.mCompositions = [];            
        end
        
        function build_mass_2_composition( this, mapfile )
            if nargin < 2, mapfile = []; end
            if this.mPermethylated
                maxMass = 3200;
            else
                maxMass = 2200;
            end
            
            n = prod( CMass2Composition.mMonoMaxNums );
            masses = zeros( n, 1 );
            compositions = zeros( n, CMonosaccharideSet.mNumberMonosaccharideClasses );
            idx = 0;
            for n1 = 0 : CMass2Composition.mMonoMaxNums(1) % 'Xyl'
                for n2 = 0 : CMass2Composition.mMonoMaxNums(2) % 'Fuc'
                    for n3 = 0 : CMass2Composition.mMonoMaxNums(3) % 'Hex'
                        for n4 = 0 : CMass2Composition.mMonoMaxNums(4) % 'HexA'
                            for n5 = 0 : CMass2Composition.mMonoMaxNums(5) % 'HexNAc'
                                for n6 = 0 : CMass2Composition.mMonoMaxNums(6) % 'Kdo'
                                    for n7 = 0 : CMass2Composition.mMonoMaxNums(7) % 'Neu5Ac'
                                        for n8 = 0 : CMass2Composition.mMonoMaxNums(8) % 'Neu5Gc'
                                            totalNum = n1 + n2 + n3 + n4 + n5 + n6 + n7 + n8;
                                            totalMass = [n1, n2, n3, n4, n5, n6, n7, n8] * this.mMonoMasses';
                                            
                                            if totalNum < 3
                                                continue;
                                            end
                                            totalMass = totalMass - (totalNum - 1) * CMass.H2O + CMass.Proton;
                                            if this.mPermethylated
                                                totalMass = totalMass - 2 * (totalNum - 1) * CMass.CH2;
                                            end
                                            if totalMass > maxMass
                                                break;
                                            end
                                            
                                            idx = idx + 1;
                                            if mod(idx, 1000) == 0
                                                fprintf('.'); 
                                                if mod(idx, 60000) == 0
                                                    fprintf(' %d\n', idx );
                                                end
                                            end
                                            masses(idx) = totalMass;
                                            compositions(idx, :) = [n1, n2, n3, n4, n5, n6, n7, n8];
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            fprintf('\n');
            
            this.mMasses = masses(1:idx);
            this.mCompositions = compositions(1:idx,:);
            [~, idxes] = sort( this.mMasses );
            this.mMasses = this.mMasses(idxes);
            this.mCompositions = this.mCompositions(idxes, :);
            
            if ~isempty(mapfile)
                save( mapfile, 'this' );
            end
        end
        
        function result = composition_to_mass( this, composition )
            n = sum( composition );
            result = this.mMonoMasses * composition' - (n-1) * CMass.H2O - 2 * (n-1) * CMass.CH2 + CMass.Proton;
        end
        
        function spectra = correct_spectrum( this, spectrum )
        % Assume spectrum.mPeaks(end) is the precursor. 
        % DO NOT add complementary peaks prior calling this function.
        % Complementary peaks are added for each possible composition in the corresponding spectrum.
            if spectrum.mProtonated == 0
                spectrum.protonate();
            end
            massAccuRange = this.mMassAccuracyPPM * 0.000001;
            
            % this.mMasses do not include the mass of the reducing-end modification
            this.set_reducing_end_modification( spectrum.mReducingEndModification );
            precursorMZ = spectrum.mPrecursor - this.mREMMassCompensation;
            if str_startswith( this.mReducingEndModification, 'REM_M3_Bion_' )
                precursorMZ = precursorMZ + CMass.H2O;
                if this.mPermethylated
                    precursorMZ = precursorMZ + CMass.CH2;
                end
            end
            this.mPrecursorMassTolerance = spectrum.mPrecursor * massAccuRange;
            
            idxes = find( abs(precursorMZ - this.mMasses) < this.mPrecursorMassTolerance );
            spectra = CSpectrum.empty(0, length(idxes));
            for k = 1 : length( idxes )
                newS = CSpectrum();
                spectra(k) = newS;

                newS.mPrecursor = this.mMasses( idxes(k) ) + this.mREMMassCompensation;
                if str_startswith( this.mReducingEndModification, 'REM_M3_Bion_' )
                    newS.mPrecursor = newS.mPrecursor - CMass.H2O;
                    if this.mPermethylated
                        newS.mPrecursor = newS.mPrecursor - CMass.CH2;
                    end
                end
                newS.mComposition = this.mCompositions(idxes(k), :);
                
                newS.mExperimentMethod = spectrum.mExperimentMethod;
                newS.mMetal = spectrum.mMetal;
                newS.mNLinked = spectrum.mNLinked;
                newS.mPermethylated = spectrum.mPermethylated;
                newS.mProtonated = spectrum.mProtonated;
                newS.mMassAccuracy = spectrum.mMassAccuracy;
                newS.mReducingEndModification = spectrum.mReducingEndModification;
                newS.comment = spectrum.comment;
                newS.filename = spectrum.filename;
                
                totalN = sum(newS.mComposition);
                v = zeros(1, totalN);
                pos = 1;
                for n = 1 : length(newS.mComposition)
                    for m = 1 : newS.mComposition(n)
                        v(pos) = n;
                        pos = pos + 1;
                    end
                end
                for n = 1 : totalN - 1
                    nreComposition = unique(nchoosek( v, n), 'row');
                    for m = 1 : size( nreComposition, 1 )
                        massB = sum( this.mMonoMasses( nreComposition(m,:) ) ) - n * CMass.H2O + CMass.Proton;
                        if this.mPermethylated
                            massB = massB - (2*n - 1) * CMass.CH2;
                        end
                        massC = massB + CMass.H2O;
                        massY = newS.mPrecursor - massB + CMass.Proton;
                        massZ = newS.mPrecursor - massC + CMass.Proton;
                    
                        peaks = newS.find_peaks(massB, 0.0001);
                        if isempty( peaks )
                            peaks = spectrum.find_peaks(massB, massB * massAccuRange);
                            if ~isempty( peaks )
                                bPeak = CPeak;
                                bPeak.mSpectrum = newS;
                                bPeak.mMass = massB;
                                bPeak.mOriPeaks = peaks;
                                bPeak.mIntensity = sum( [peaks.mIntensity] );
                                bPeak.mMassRange = [massB-0.0001, massB+0.0001];
                                bPeak.mType = CPeakType.B;
                                newS.mPeaks(end+1) = bPeak;
                            else % If Y-ion exists, add complementary ion.
                                peaks = spectrum.find_peaks(massY, massY * massAccuRange);
                                if ~isempty( peaks )
                                    bPeak = CPeak;
                                    bPeak.mSpectrum = newS;
                                    bPeak.mMass = massB;
                                    bPeak.mIsComplement = 1;
                                    bPeak.mOriPeaks = peaks;
                                    bPeak.mIntensity = sum( [peaks.mIntensity] );
                                    bPeak.mMassRange = [massB-0.0001, massB+0.0001];
                                    bPeak.mType = CPeakType.B;
                                    newS.mPeaks(end+1) = bPeak;
                                end
                            end
                        else % The peak already exists.
                            for h = 1 : length(peaks)
                                peaks(h).mType = CPeakType.combineType( peaks(h).mType, CPeakType.B );
                            end
                        end
                        if this.mCheckMinus2H % Do this only for B- and C- ions.
                            massB = massB - CMass.H2;
                            peaks = newS.find_peaks(massB, 0.0001);
                            if isempty( peaks )
                                peaks = spectrum.find_peaks(massB, massB * massAccuRange);
                                if ~isempty( peaks )
                                    bPeak = CPeak;
                                    bPeak.mSpectrum = newS;
                                    bPeak.mMass = massB;
                                    bPeak.mOriPeaks = peaks;
                                    bPeak.mIntensity = sum( [peaks.mIntensity] );
                                    bPeak.mMassRange = [massB-0.0001, massB+0.0001];
                                    bPeak.mType = CPeakType.combineType(CPeakType.B, CPeakType.Minus2H);
                                    newS.mPeaks(end+1) = bPeak;
                                end
                            else
                                for h = 1 : length(peaks)
                                    peaks(h).mType = CPeakType.combineType( peaks(h).mType, CPeakType.combineType(CPeakType.B, CPeakType.Minus2H) );
                                end
                            end
                        end

                        peaks = newS.find_peaks(massC, 0.0001);
                        if isempty( peaks )
                            peaks = spectrum.find_peaks(massC, massC * massAccuRange);
                            if ~isempty( peaks )
                                cPeak = CPeak;
                                cPeak.mSpectrum = newS;
                                cPeak.mMass = massC;
                                cPeak.mOriPeaks = peaks;
                                cPeak.mIntensity = sum( [peaks.mIntensity] );
                                cPeak.mMassRange = [massC-0.0001, massC+0.0001];
                                cPeak.mType = CPeakType.C;
                                newS.mPeaks(end+1) = cPeak;
                            else % If Z-ion exists, add complementary ion.
                                peaks = spectrum.find_peaks(massZ, massZ * massAccuRange);
                                if ~isempty( peaks )
                                    cPeak = CPeak;
                                    cPeak.mSpectrum = newS;
                                    cPeak.mMass = massC;
                                    cPeak.mIsComplement = 1;
                                    cPeak.mOriPeaks = peaks;
                                    cPeak.mIntensity = sum( [peaks.mIntensity] );
                                    cPeak.mMassRange = [massC-0.0001, massC+0.0001];
                                    cPeak.mType = CPeakType.C;
                                    newS.mPeaks(end+1) = cPeak;
                                end
                            end
                        end
                        if this.mCheckMinus2H % Do this only for B- and C- ions.
                            massC = massC - CMass.H2;
                            peaks = newS.find_peaks(massC, 0.0001);
                            if isempty( peaks )
                                peaks = spectrum.find_peaks(massC, massC * massAccuRange);
                                if ~isempty( peaks )
                                    cPeak = CPeak;
                                    cPeak.mSpectrum = newS;
                                    cPeak.mMass = massC;
                                    cPeak.mOriPeaks = peaks;
                                    cPeak.mIntensity = sum( [peaks.mIntensity] );
                                    cPeak.mMassRange = [massC-0.0001, massC+0.0001];
                                    cPeak.mType = CPeakType.combineType(CPeakType.C, CPeakType.Minus2H);
                                    newS.mPeaks(end+1) = cPeak;
                                end
                            else
                                for h = 1 : length(peaks)
                                    peaks(h).mType = CPeakType.combineType( peaks(h).mType, CPeakType.combineType(CPeakType.C, CPeakType.Minus2H) );
                                end
                            end
                        end
                    end % for m = 1 : size( nreComposition, 1 )
                end % for n = 1 : totalN - 1
                
                % add complementary peaks
                n = length(newS.mPeaks);
                for m = 1 : n
                    aPeak = newS.mPeaks(m);
                    if CPeakType.isB(aPeak.mType) || CPeakType.isC(aPeak.mType)
                        if aPeak.mComplement ~= 0
                            continue;
                        end
                        for c = m+1 : n
                            if abs(aPeak.mMass + newS.mPeaks(c).mMass - newS.mPrecursor - CMass.Proton) < 0.000001
                                aPeak.mComplement = c;
                                newS.mPeaks(c).mComplement = m;
                            end
                        end
                    else
                        tempMass = newS.mPrecursor - aPeak.mMass + CMass.Proton;
                        if tempMass < this.mMonoMasses(1) - CMass.CH2 - CMass.H2O % not too small.
                            continue;
                        end
                        peaks = newS.find_peaks(tempMass, 0.0001);
                        if isempty( peaks )
                            newPeak = CPeak;
                            newPeak.mSpectrum = newS;
                            newPeak.mOriPeaks = aPeak;
                            newPeak.mMass = tempMass;
                            newPeak.mMassRange = [tempMass-0.0001, tempMass+0.0001];
                            newPeak.mIsComplement = 1;
                            newPeak.mComplement = -m; % index of aPeak
                            if CPeakType.isY( aPeak.mType )
                                newPeak.mType = CPeakType.B;
                            end
                            if CPeakType.isZ( aPeak.mType )
                                newPeak.mType = CPeakType.combineType( newPeak.mType, CPeakType.C );
                            end
                            newS.mPeaks(end+1) = newPeak;
                            aPeak.mComplement = newS.num;
                        end
                    end
                end
                newS.sort_peaks();
                
                % Add the precuror
                newPeak = CPeak;
                newPeak.mSpectrum = newS;
                newPeak.mMass = newS.mPrecursor;
                newPeak.mType = CPeakType.T;
                newPeak.mMassRange = [newPeak.mMass - 0.0001, newPeak.mMass + 0.0001];
                newPeak.mOriPeaks = spectrum.mPeaks(end);
                newS.mPeaks(end+1) = newPeak;
            end
        end
        
        function this = load( this, datafile )
            temp = load(datafile);
            this.mMonoMasses = temp.this.mMonoMasses;
            this.mCheckMinus2H = temp.this.mCheckMinus2H;
            this.mMassAccuracyPPM = temp.this.mMassAccuracyPPM;
            this.mPrecursorMassTolerance = temp.this.mPrecursorMassTolerance;
            this.mMasses = temp.this.mMasses;
            this.mCompositions = temp.this.mCompositions;
        end
        
        function result = get_permethylation( this )
            result = this.mPermethylated;
        end
        
        function set_permethylation( this, newPermethylated )
            if isempty(this.mPermethylated) || (this.mPermethylated ~= newPermethylated)
                this.mPermethylated = newPermethylated;
                for k = 1 : CMonosaccharideSet.cNumberMonosaccharideClasses
                    this.mMonoMasses(k) = CMonosaccharideSet.get_mass( CMonosaccharideSet.cMonoClasses{k}, newPermethylated );
                end
            end
        end
        
        function result = get_reducing_end_modification( this )
            result = this.mReducingEndModification;
        end
        
        function set_reducing_end_modification( this, REM )
            if strcmp(this.mReducingEndModification, REM) ~= 1
                this.mReducingEndModification = REM;
                this.mREMMassCompensation = CMass.get_mass_compensation( REM, this.mPermethylated );
            end
        end
        
        function result = get_REM_mass_compensation( this )
            result = this.mREMMassCompensation;
        end
    end
end
