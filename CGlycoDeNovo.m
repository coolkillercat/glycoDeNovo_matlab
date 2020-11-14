classdef CGlycoDeNovo < handle % By Pengyu Hong @ Brandeis University
    properties (Constant)
        cMonosaccharideClasses = { 'Xyl', 'Fuc', 'Hex', 'HexA', 'HexNAc', 'Kdo', 'Neu5Ac', 'Neu5Gc'};
        cCleavages = CCrossRingCleavage.cCRC;
        cLegalGlycosidicBonds = ...
              [1 1 1 1 1 1 1 1; ...
               1 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1; ...
               0 0 1 1 1 1 1 1];
    end
    
    properties
        mCheckMinus2H = 0;
        mCheckMinusH = 0;
        mCheckGap = 0; % if mCheckGap, allows one gap.
        mMassAccuracyPPM = 5; % ppm = 1000000*mass_error/exact_mass; mass_error = mMassAccuracyPPM*exact_mass/1000000;
        mMassAccuracyDalton = 0.005;
        mMaxBranchingNum = 2; % default bi-branching
        mIonMass = 0;
        mIonMetal = [];
        
        mReducingEndModification = '';
        
        mIntensityThreshold = 0;
        mNLinked = 0;
        mPermethylated = 0;
        mUseComplementOnly = 0
        
        mPeaks = CPeak.empty(1,0);
        mPeakWeights = [];
        
        mPossibleMonosaccharideClasses = [];
        mPossibleMonosaccharideClassIDs = [];
        mReducingEndCompensation = 0;
        mFinalPeakCompensation = 0;
        mCompositionConstraint = zeros(1, CMonosaccharideSet.cNumberMonosaccharideClasses);
        
        mTopologySuperSets = CTopologySuperSet.empty(0,0);
    end
    
    properties (SetAccess = private)
        mDelta = [];
        mDelta2 = [];

        mNumPeaks = 0;
        mCurrentPeakIdx = 0;
        mCurrentTargetMass = 0;
        mCurrentTargetMassLow = 0;
        mCurrentTargetMassHigh = 0;
        mCurrentTPSuperSetSize = 0;
        mCurrentBranches = {[], [], [], []};
        mLeafPeak = 0;
        mFinalPeak = 0;
        mPeakMasses = [];
        mLinkageScores = [];
        mLinkageMatchedPeakIdx = {};
        mMassAccuracyPP = 0.000005;
        
        mTryCIon = 0;
        mCurrentTopologySuperSetB = CTopologySuperSet.empty(0,0);
        mCurrentTopologySuperSetC = CTopologySuperSet.empty(0,0);
        mCurrentTopologySuperSetB1 = CTopologySuperSet.empty(0,0);
        mCurrentTopologySuperSetC1 = CTopologySuperSet.empty(0,0);
        mCurrentTopologySuperSetB2 = CTopologySuperSet.empty(0,0);
        mCurrentTopologySuperSetC2 = CTopologySuperSet.empty(0,0);
        
        mPermMatrixes = cell(3,4);
        
        %mMeetMonoCountThreshold = Inf;
        
        mGlycan = [];
        mGlycanUnits = []; % mGlycanUnits(k) is the k-th unit
        mGlycanChildrenMap = []; % mGlycanChildrenMap(k) is the children of the k-th unit
        mGlycanUnitLinkedInMap = [];  % mGlycanUnitLinkedInMap(k) contains all possible linkedIn configurations to the k-th unit
        mGlycanUnitLinkedInMapIdx = []; % mGlycanUnitLinkedInMap(mGlycanUnitLinkedInMapIdx(k), :) is the linkedIn configuration of the k-th unit 
        mGlycanUnitCleavageMassMap = []; % mGlycanUnitCleavageMassMap{k}{m} = the masses of ions by cleaving the k-th unit with the m-th linkedIn configuration
        mGlycanUnitCleavageIonMap = []; % mGlycanUnitCleavageMassMap{k}{m} = the ions of ions by cleaving the k-th unit with the m-th linkedIn configuration
        
        mPotentialSpectrumSet = {};
        mSpectrumScoreMatrix = [];
        mCompleteSpectrum = [];
        mTotalNumConfiguration = 0;
    end

    methods (Static)
        function flag = isLegalGlycosidicBond( ion_or_mono, toMonoClassID )
        % find out if "ion_or_mono -> mono" is a legal glycosidic bond.
            if isa( ion_or_mono, 'double' ) && ( ion_or_mono >= 1 ) && ( ion_or_mono <= 8 )
                fromMonoClassID = ion_or_mono;
            elseif isa( ion_or_mono, 'CIon' )
                if isempty( ion_or_mono ) || isempty( ion_or_mono.stem )
                    fromMonoClassID = 0;
                else
                    fromMonoClassID = ion_or_mono.stem(1).mClassID;
                end
            elseif isa( ion_or_mono, 'CMonosaccharide' )
                fromMonoClassID = ion_or_mono.mClassID;
            else
                throw( MException('CGlycoDeNovo:isLegalGlycosidicBond', 'ion_or_mono is not a recognizable input argument.') );
            end
            
            flag = 0;
            if fromMonoClassID == 0
                if (toMonoClassID > 0) && (toMonoClassID < 8)
                    flag = 1;
                end
            else
                flag = CGlycoDeNovo.cLegalGlycosidicBonds(fromMonoClassID, toMonoClassID);
            end
        end
    end
        
    methods
        function obj = CGlycoDeNovo( massAccuracy, possibleMonoClasses, minus2H, allowGap )
            if nargin < 1 || isempty( massAccuracy )
                obj.mMassAccuracyPPM = 5; % ppm
            else
                obj.mMassAccuracyPPM = massAccuracy;
            end
            obj.mMassAccuracyPP = obj.mMassAccuracyPPM / 1000000;
            if nargin < 2
                obj.mPossibleMonosaccharideClasses = [];
                obj.mPossibleMonosaccharideClassIDs = [];
            else
                if iscell( possibleMonoClasses )
                    obj.mPossibleMonosaccharideClasses = possibleMonoClasses;
                    obj.mPossibleMonosaccharideClassIDs = CMonosaccharideSet.find_monosaccharide_classID(possibleMonoClasses);
                else
                    obj.mPossibleMonosaccharideClasses = CMonosaccharideSet.find_className_by_classID(possibleMonoClasses);
                    obj.mPossibleMonosaccharideClassIDs = possibleMonoClasses;
                end
            end
            if nargin < 3 || isempty( minus2H )
                obj.mCheckMinus2H = 0;
            else
                obj.mCheckMinus2H = minus2H;
            end
            if nargin < 4 || isempty( allowGap )
                obj.mCheckGap = 0;
            else
                obj.mCheckGap = allowGap;
            end
            
            obj.mCompositionConstraint = ones(1,8) * 1000;
            
            % prepare mPermMatrixes, which will be used in scoring linkages
            for b = 1 : 3
                for numV = 2 : 4
                    if numV < b, continue; end
                    temp = perms(1:numV);
                    obj.mPermMatrixes{b, numV} = unique( temp(:, 1:b), 'rows' );
                end
            end
        end

        function init( obj, compositionConstraint )
            %obj.mReducingEndCompensation = CMass.get_mass_compensation( obj.mReducingEndModification, obj.mPermethylated );
            %obj.mFinalPeakCompensation = obj.mReducingEndCompensation + obj.mPermethylated * CMass.CH2 + CMass.H2O;
            
            if nargin == 2
                assert( ~any(compositionConstraint < 0) )
                obj.mCompositionConstraint = compositionConstraint;
            end
            
            perMassLoss = 2 * obj.mPermethylated * CMass.CH2;  % If no the root unit, lose two CH2s when permethylated

            linkageMassLoss = CMass.H2O + perMassLoss;
            numClass = length( CGlycoDeNovo.cMonosaccharideClasses );
            obj.mDelta = [];
            for k = 1 : numClass
                mono = CMonosaccharide( CGlycoDeNovo.cMonosaccharideClasses{k}, obj.mPermethylated );
                mono.update();

                % add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a B-ion
                obj.mDelta.B2B.mass(k) = mono.mMass - linkageMassLoss;
                obj.mDelta.B2B.unit(k) = mono; %mono.copy; 
                
                % add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a C-ion
                obj.mDelta.B2C.mass(k) = mono.mMass - perMassLoss;
                obj.mDelta.B2C.unit(k) = mono; % mono.copy; 
                
                % add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a B-ion
                obj.mDelta.C2B.mass(k) = mono.mMass - CMass.H2O - linkageMassLoss;
                obj.mDelta.C2B.unit(k) = mono; % mono.copy; 
                
                % add 1 CGlycoDeNovo.cMonosaccharideClasses(k) to form a C-ion
                obj.mDelta.C2C.mass(k) = mono.mMass - CMass.H2O - perMassLoss;
                obj.mDelta.C2C.unit(k) = mono; % mono.copy; 
            end
            
            obj.mDelta2 = [];
            idx = 1;
            for k = 1 : numClass
                mono1 = CMonosaccharide( CGlycoDeNovo.cMonosaccharideClasses{k}, obj.mPermethylated );
                for m = 1 : numClass
                    mono2 = CMonosaccharide( CGlycoDeNovo.cMonosaccharideClasses{m}, obj.mPermethylated );

                    if CGlycoDeNovo.isLegalGlycosidicBond( mono1, mono2.mClassID ) == 0
                        continue;
                    end
                    
                    % B-ion -> mono1 -> mono2 ==> B-ion
                    obj.mDelta2.B2B.mass(idx) = mono1.mMass + mono2.mMass - linkageMassLoss*2;
                    obj.mDelta2.B2B.unit(idx, 1) = mono1.copy;
                    obj.mDelta2.B2B.unit(idx, 2) = mono2.copy;
                    
                    % C-ion -> mono1 -> mono2 ==> B-ion
                    obj.mDelta2.B2C.mass(idx) = mono1.mMass + mono2.mMass - linkageMassLoss + CMass.H2O;
                    obj.mDelta2.B2C.unit(idx, 1) = mono1.copy;
                    obj.mDelta2.B2C.unit(idx, 2) = mono2.copy;
                    
                    % B-ion -> mono1 -> mono2 ==> C-ion
                    obj.mDelta2.C2B.mass(idx) = mono1.mMass + mono2.mMass - CMass.H2O - linkageMassLoss * 2;
                    obj.mDelta2.C2B.unit(idx, 1) = mono1.copy;
                    obj.mDelta2.C2B.unit(idx, 2) = mono2.copy;
                    
                    % C-ion -> mono1 -> mono2 ==> C-ion
                    obj.mDelta2.C2C.mass(idx) = mono1.mMass + mono2.mMass - linkageMassLoss * 2;
                    obj.mDelta2.C2C.unit(idx, 1) = mono1.copy;
                    obj.mDelta2.C2C.unit(idx, 2) = mono2.copy;
                    
                    idx = idx + 1;
                end
            end
        end

        function pos = acceptable_monosaccharide( obj, monoClass )
            if isempty( obj.mPossibleMonosaccharideClasses )
                pos = 1;
            else
                pos = sum( strcmp( obj.mPossibleMonosaccharideClasses, monoClass ) ) > 0;
            end
        end

        function interpret_peaks( obj, spec )
        % spec.peaks are arranaged from lightest to heaviest.
        % spec.peaks may need to be pre-processed to include complement ions.
            if isempty( spec ) || isempty( spec.mPeaks )
                return; 
            end
            
            obj.mMassAccuracyDalton = obj.mMassAccuracyPP * spec.mPrecursor;
            if spec.mProtonated
                obj.mIonMetal = 'Proton'; 
                obj.mIonMass = CMass.Proton;
            else
                obj.mIonMetal = spec.mMetal;
                switch spec.mMetal
                    case 'H'
                        obj.mIonMetal = 'Proton';
                        obj.mIonMass = CMass.Proton;
                    case 'Proton'
                        obj.mIonMass = CMass.Proton;
                    case 'Lithium'
                        obj.mIonMass = CMass.Lithium - CMass.Electron;
                    case 'Na'
                        obj.mIonMass = CMass.Sodium - CMass.Electron;
                    case 'Sodium'
                        obj.mIonMass = CMass.Sodium - CMass.Electron;
                    otherwise
                        ME = MException('CGlycoDeNovo:CGlycoDeNovo', 'ionMetal %s not found',ionMetal);
                        throw(ME);
                end
            end
            obj.mPermethylated = spec.mPermethylated;
            obj.mNLinked = spec.mNLinked;
            obj.mReducingEndModification = spec.mReducingEndModification;
            obj.init();
        
            spec.clear_inferred();
            obj.mPeaks = spec.mPeaks;
            obj.mPeakMasses = [obj.mPeaks.mMass];
            obj.mNumPeaks = length( spec.mPeaks );
            obj.mCurrentPeakIdx = 0;
            % obj.mMeetMonoCountThreshold = length( obj.mCompositionConstraint );
            
            obj.mTopologySuperSets = CTopologySuperSet.empty(0,0); % add an empty super set
            
            endFlank = ' ---';
            fprintf( '--- Peak ' );
            obj.mFinalPeak = 0;
            for k = 1 : obj.mNumPeaks
                if obj.mPeaks(k).mIntensity < obj.mIntensityThreshold, continue; end
                
                % Per request by Cheng Lin 2017/07. Test using
                % complementary peaks only. Did not work well.
                if (k ~= obj.mNumPeaks) && obj.mUseComplementOnly && ...
                   ~obj.mPeaks(k).mIsComplement
                   % (~obj.mPeaks(k).mIsComplement && (obj.mPeaks(k).mComplement <= 0) )
                    continue;
                end
                
                obj.mCurrentTargetMass = obj.mPeaks(k).mMass;
                obj.mCurrentTargetMassLow = obj.mPeaks(k).mMassRange(1);
                obj.mCurrentTargetMassHigh = obj.mPeaks(k).mMassRange(2);
                obj.mCurrentTopologySuperSetC = [];
                obj.mCurrentTopologySuperSetB = [];
                obj.mCurrentTopologySuperSetB1 = [];
                obj.mCurrentTopologySuperSetC1 = [];
                obj.mCurrentTopologySuperSetB2 = [];
                obj.mCurrentTopologySuperSetC2 = [];
                
                if obj.mPermethylated
                    if obj.mCurrentTargetMass < 175, continue; end % too light
                else
                    if obj.mCurrentTargetMass < 131, continue; end % too light
                end
                
                obj.mCurrentPeakIdx = k;
                obj.mCurrentTPSuperSetSize = length( obj.mTopologySuperSets );
                obj.mCurrentBranches = {[], [], [], []};

                msg = [num2str(k), '/', num2str(obj.mNumPeaks), ' (', ...
                       num2str(length(obj.mTopologySuperSets)), ')', endFlank];
                fprintf( '%s', msg );

                if obj.mCurrentPeakIdx == obj.mNumPeaks
                    obj.mTryCIon = 0;
                    obj.mFinalPeak = 1;
                else
                    obj.mTryCIon = 1;
                    
                    bionMass = obj.mCurrentTargetMass - CMass.H2O;
                    for m = obj.mCurrentTPSuperSetSize : -1 : 1
                        temp = (obj.mTopologySuperSets(m).mMassPros(2) + obj.mTopologySuperSets(m).mMassPros(3))/2;
                        if bionMass > temp + 0.01
                            break;
                        elseif bionMass >= temp - 0.01
                            obj.mTopologySuperSets(m).add_peak( obj.mCurrentPeakIdx, 2 );
                            obj.mPeaks(k).mInferredSuperSet = obj.mTopologySuperSets(m);
                            obj.mTryCIon = 0;
                            break;
                        end
                    end
                end
                
                % obj.interpret_a_peak_recursion( 1, 0 );
                obj.interpret_a_peak();

                if k == obj.mNumPeaks && isempty( obj.mCurrentTopologySuperSetB ) && isempty( obj.mCurrentTopologySuperSetC )
                    obj.append_NLinkedRoot();
                end
                
                obj.add_currentTSS_to_pool();
                
                if k < obj.mNumPeaks
                    for m = 1 : length( msg )
                        fprintf( '\b' );
                    end
                else
                    fprintf( '\n' );
                end
            end
            
            disp( 'CGlycoDeNovo::reconstruct_topology done!' );
        end
        
        function interpret_a_peak_recursion( obj, checkBranchingIdx, startTPSuperSetIdx ) % This is significantly slower in Matlab. Maybe Java/C implementation can be faster.
            if ( startTPSuperSetIdx == 0 ) % leaf cleavages, happens only @ checkBranchingIdx = 1
               if ( obj.mCurrentTargetMass < 438 )  % should not be heavier than the heaviest monosaccharide
                  obj.test_topologyset( 1 );
               end
               startTPSuperSetIdx = 1;
            end
            
            for k = startTPSuperSetIdx : obj.mCurrentTPSuperSetSize
                obj.mCurrentBranches{checkBranchingIdx} = obj.mTopologySuperSets(k);
                result = obj.test_topologyset( checkBranchingIdx );
                if result == -1 % no need to try further because too heavy.
                    % obj.mCurrentBranches{checkBranchingIdx} = [];
                    break;
                end
                if checkBranchingIdx < obj.mMaxBranchingNum
                    obj.interpret_a_peak_recursion( checkBranchingIdx+1, k ); % add more branches, and try.
                end
            end
            obj.mCurrentBranches{checkBranchingIdx} = [];
        end
        
        function interpret_a_peak( obj )
            
            if ( obj.mCurrentTargetMass < 438 )  % Interpret as a monosaccharide. The mass should not be heavier than the heaviest monosaccharide
                if obj.mPermethylated
                    massCompensationLow = CMass.CH2 + obj.mIonMass; % CMass.CH2 was pre-subtracted because obj.mDelta is prepared for non-end units.
                else
                    massCompensationLow = obj.mIonMass; % CMass.CH2 was pre-subtracted because obj.mDelta is prepared for non-end units.
                end
                if obj.mFinalPeak
                    massCompensationLow = massCompensationLow + obj.mFinalPeakCompensation;
                end
                massCompensationHigh = massCompensationLow * (1 + obj.mMassAccuracyPP);
                massCompensationLow = massCompensationLow * (1 - obj.mMassAccuracyPP);
                obj.mLeafPeak = 1;
                result = obj.test_topologyset( [massCompensationLow, massCompensationHigh] );
                if result == -1 % no need to try further because too heavy.
                    return;
                end
            end
            
            branchMasses = zeros(4,2);
            massD = obj.mIonMass + obj.mPermethylated * CMass.CH2; % each branch causes a CH2 loss to the joint monosaccharide when permethylated && obj.mDelta only considers linear structure (i.e., one branch).;
            obj.mLeafPeak = 0;
            for k = 1 : obj.mCurrentTPSuperSetSize  % 1st branch
                obj.mCurrentBranches{1} = obj.mTopologySuperSets(k);
                if obj.mFinalPeak
                    branchMasses(1,:) = obj.mCurrentBranches{1}.mMassPros(2:3) + obj.mFinalPeakCompensation * [1 - obj.mMassAccuracyPP, 1 + obj.mMassAccuracyPP];
                else
                    branchMasses(1,:) = obj.mCurrentBranches{1}.mMassPros(2:3);
                end
                result = obj.test_topologyset( branchMasses(1,:) );
                if result == -1 % no need to try further because too heavy.
                    break;
                end
                
                for kk = k : obj.mCurrentTPSuperSetSize  % 2nd branch
                    obj.mCurrentBranches{2} = obj.mTopologySuperSets(kk);
                    branchMasses(2,:) = obj.mCurrentBranches{2}.mMassPros(2:3) - massD;
                    result = obj.test_topologyset( sum(branchMasses) );
                    if result == -1 % no need to try further because too heavy.
                        break;
                    end
                    
                    % S: 2018/12/29. 2-branch-with-gap
                    if (obj.mMaxBranchingNum < 3) && (obj.mCheckGap == 0) % maximum number of branches = 2
                        continue; 
                    end
                    % E: 2018/12/29. 2-branch-with-gap
                    
                    for kkk = kk : obj.mCurrentTPSuperSetSize  % 3rd branch
                        obj.mCurrentBranches{3} = obj.mTopologySuperSets(kkk);
                        branchMasses(3,:) = obj.mCurrentBranches{3}.mMassPros(2:3) - massD * [1 + obj.mMassAccuracyPP, 1 - obj.mMassAccuracyPP];
                        
                        % S: 2018/12/29. 2-branch-with-gap
                        if (obj.mCheckGap > 0) && (obj.mMaxBranchingNum < 3)
                            result = obj.test_topologyset_2BranchWithGap( sum(branchMasses) );
                        else
                            result = obj.test_topologyset( sum(branchMasses) );
                        end
                        % E: 2018/12/29. 2-branch-with-gap
                        
                        if result == -1 % no need to try further because too heavy.
                            break;
                        end
                        
                        if obj.mMaxBranchingNum < 4 % maximum number of branches = 3
                            continue;
                        end
                        
                        for kkkk = kkk : obj.mCurrentTPSuperSetSize  % 4th branch
                            obj.mCurrentBranches{4} = obj.mTopologySuperSets(kkkk);
                            branchMasses(4,:) = obj.mCurrentBranches{4}.mMassPros(2:3) - massD * [1 - obj.mMassAccuracyPP, 1 + obj.mMassAccuracyPP];;
                            result = obj.test_topologyset( sum(branchMasses) );
                            if result == -1 % no need to try further because too heavy.
                                break;
                            end
                        end
                        obj.mCurrentBranches{4} = [];
                        branchMasses(4,:) = 0;
                    end
                    obj.mCurrentBranches{3} = [];
                    branchMasses(3,:) = 0;
                end
                obj.mCurrentBranches{2} = [];
                branchMasses(2,:) = 0;
            end
            obj.mCurrentBranches{1} = [];
        end
        
        function set_reducing_end_modification( this, REM )
            this.mReducingEndModification = REM;
            % this.mReducingEndCompensation = CGlycan.reducing_end_mass_compensation( this.mReducingEndModification, this.mPermethylated );
            this.mReducingEndCompensation = CMass.get_mass_compensation( this.mReducingEndModification, this.mPermethylated );
            if str_startswith( REM, 'REM_M3_Bion_' )
                this.mFinalPeakCompensation = this.mReducingEndCompensation; % CH2 and H2O already included.
            elseif ~str_startswith( REM, 'REM_M3_' )
                this.mFinalPeakCompensation = this.mReducingEndCompensation + this.mPermethylated * CMass.CH2 + CMass.H2O;
            end
        end
        
        function result = test_topologyset( obj, massBound )
            result = 0;
            massCompensationLow = massBound(1);
            massCompensationHigh = massBound(2);
            
            % Check if the branches together are already too heavy or still too light
            temp = obj.mCurrentTargetMass - massCompensationLow;
            if obj.mPermethylated
                if temp < 160
                    result = -1; return;
                end
                if obj.mCheckGap
                    if temp > 860 % max( obj.mDelta2.B2C )
                        return;
                    end
                elseif temp > 420 % max( obj.mDelta.B2C )
                    return;
                end
            else
                if temp < 105
                    result = -1; return;
                end
                if obj.mCheckGap
                    if temp > 660 % max( obj.mDelta2.B2C )
                        return;
                    end
                elseif temp > 335 % max( obj.mDelta.B2C )
                    return;
                end
            end
            
            % (2) try to extend to B-ions by adding one monosaccharide 
            % @@@@ Here. Need to deal with [massCompensationLow, massCompensationHigh]
            theoreticalMass = obj.mDelta.B2B.mass + [massCompensationLow; massCompensationHigh];
            %theoreticalMass = [obj.mDelta.B2B.mass * (1 - obj.mMassAccuracyPP) + massCompensationLow; ...
            %                   obj.mDelta.B2B.mass * (1 + obj.mMassAccuracyPP) + massCompensationHigh];
            idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                        ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
            flagMinusH = zeros(1, length(idx));
            if ~obj.mFinalPeak
                if obj.mCheckMinus2H
                    temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                 ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                    idx = [idx, temp];
                    flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                end
                if obj.mCheckMinusH
                    temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                 ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                    idx = [idx, temp];
                    flagMinusH = [flagMinusH, ones(1, length(temp))];
                end
            end
            
            for a = 1 : length(idx)
                newUnit = obj.mDelta.B2B.unit( idx(a) );
%                if ( obj.acceptable_monosaccharide( newUnit.mClass ) == 0 ) || ...
                if ~any( obj.mPossibleMonosaccharideClassIDs == newUnit.mClassID ) || ...
                   ( newUnit.mClassID == 2 && obj.mLeafPeak == 0 ) || ...
                   ( obj.mCompositionConstraint( newUnit.mClassID ) < 1 )
                    continue;
                end
                
                newSet = CTopologySet;
                temp = flagMinusH(a)*CMass.H;
                newSet.mMasses = [obj.mCurrentTargetMass, ...
                                  max(theoreticalMass(1, idx(a)), obj.mCurrentTargetMassLow + temp), ...
                                  min(theoreticalMass(2, idx(a)), obj.mCurrentTargetMassHigh + temp)];
                newSet.mRootMono = newUnit;
                newSet.mSources = obj.mCurrentBranches;
                newSet.mReconstructor = obj;
                newSet.mMinusH = flagMinusH(a);
                if obj.mFinalPeak
                    newSet.mType = 'T';
                    if flagMinusH(a) == 0
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T];
                    elseif flagMinusH(a) == 1
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T_H];
                    else
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T_2H];
                    end
                else
                    newSet.mType = 'B';
                    if flagMinusH(a) == 0
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B];
                    elseif flagMinusH(a) == 1
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B_H];
                    else
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B_2H];
                    end
                end
                obj.insert_into_currentTopologySuperSet( newSet );
                % S: 2018/12/29. 2-branch-with-gap
                result = 1;
                % E: 2018/12/29. 2-branch-with-gap
            end
            
            failed = isempty(idx);
            
            if ~obj.mFinalPeak && obj.mTryCIon % not the precursor && allowed C-ion, then try to form C-Ions
                theoreticalMass = obj.mDelta.B2C.mass + [massCompensationLow; massCompensationHigh];
                %theoreticalMass = [obj.mDelta.B2C.mass * (1 - obj.mMassAccuracyPP) + massCompensationLow; 
                %                   obj.mDelta.B2C.mass * (1 + obj.mMassAccuracyPP) + massCompensationHigh];
                idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                            ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
                flagMinusH = zeros(1, length(idx));
                if obj.mCheckMinus2H
                    temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                 ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                    idx = [idx, temp];
                    flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                end
                if obj.mCheckMinusH
                    temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                 ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                    idx = [idx, temp];
                    flagMinusH = [flagMinusH, ones(1, length(temp))];
                end
                
                for a = 1 : length(idx)
                    newUnit = obj.mDelta.B2C.unit(idx(a));
                    if ( obj.acceptable_monosaccharide( newUnit.mClass ) == 0 ) || ...
                       ( newUnit.mClassID == 2 && ~obj.mLeafPeak ) || ...
                       ( obj.mCompositionConstraint( newUnit.mClassID ) < 1 )
                        continue;
                    end
                    
                    newSet = CTopologySet;
                    temp = flagMinusH(a)*CMass.H;
                    
                    newSet.mMasses = [obj.mCurrentTargetMass, ...
                                      [max(theoreticalMass(1, idx(a)), obj.mCurrentTargetMassLow + temp), ...
                                       min(theoreticalMass(2, idx(a)), obj.mCurrentTargetMassHigh + temp)] - CMass.H2O];  % Although it is C-ion, we use its B-ion mass.
                    newSet.mRootMono = newUnit;
                    newSet.mSources = obj.mCurrentBranches;
                    newSet.mType = 'C';
                    newSet.mReconstructor = obj;
                    newSet.mMinusH = flagMinusH(a);
                    % newSet.mTargetPeaks = [obj.mCurrentPeakIdx; 2]; % Edited: 2019/02/15.
                    if newSet.mMinusH == 0 % Edited: 2020/06/23.
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C]; 
                    elseif newSet.mMinusH == 1
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C_H];
                    else
                        newSet.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C_2H];
                    end
                    
                    obj.insert_into_currentTopologySuperSet( newSet );
                    % S: 2018/12/29. 2-branch-with-gap
                    result = 1;
                    % E: 2018/12/29. 2-branch-with-gap
                end
                
                failed = failed && isempty(idx);
            end

            % Allow 1 gap only when no-gap fails
            if (failed && obj.mCheckGap) || (obj.mCheckGap == 2)
                if isempty( obj.mCurrentBranches{2} ) && ~isempty( obj.mCurrentBranches{1} ) % Rule: Fuc can only be a branch, in a linear substructure.
                    possible = 0;
                    for aTPS = obj.mCurrentBranches{1}.mTopologySets
                        if aTPS.mRootMono.mClassID ~= 2
                            possible = 1;
                            break;
                        end
                    end
                    if ~possible, return; end
                end
                
                % first try B-ions
                theoreticalMass = obj.mDelta2.B2B.mass + [massCompensationLow; massCompensationHigh];
                %theoreticalMass = [obj.mDelta2.B2B.mass * (1 - obj.mMassAccuracyPP) + massCompensationLow; ...
                %                   obj.mDelta2.B2B.mass * (1 + obj.mMassAccuracyPP) + massCompensationHigh];
                idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                            ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
                flagMinusH = zeros(1, length(idx));
                if ~obj.mFinalPeak
                    if obj.mCheckMinus2H
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                    end
                    if obj.mCheckMinusH
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))];
                    end
                end
                
                for a = 1 : length(idx)
                    newUnit1 = obj.mDelta2.B2B.unit(idx(a), 1);
                    if newUnit1.mClassID == 2  % Rule: Fuc can only be a branch, not in a linear substructure.
                        continue;
                    end
                    newUnit2 = obj.mDelta2.B2B.unit(idx(a), 2);
                    % if obj.acceptable_monosaccharide( newUnit1.mClass ) == 0 || ...
                    %   obj.acceptable_monosaccharide( newUnit2.mClass ) == 0
                    %    continue;
                    % end

                    temp = flagMinusH(a)*CMass.H;
                    massPros = [obj.mCurrentTargetMass, ...
                                max(theoreticalMass(1, idx(a)), obj.mCurrentTargetMassLow + temp), ...
                                min(theoreticalMass(2, idx(a)), obj.mCurrentTargetMassHigh + temp)];

                    newSet = CTopologySet;
                    newSet.mMasses = massPros;
                    newSet.mRootMono = newUnit2;
                    newSet.mMissingMono = newUnit1; % Add gap unit to the stem.
                    newSet.mSources = obj.mCurrentBranches;
                    newSet.mReconstructor = obj;
                    newSet.mMinusH = flagMinusH(a);
                    
                    if obj.mFinalPeak
                        newSet.mType = 'T';
                    else
                        newSet.mType = 'B';
                    end
                    obj.insert_into_currentTopologySuperSet( newSet );
                    % S: 2018/12/29. 2-branch-with-gap
                    result = 1;
                    % E: 2018/12/29. 2-branch-with-gap
                    
                    for b = 1 : 4
                        if isempty( obj.mCurrentBranches{b} )
                            break;
                        end
                        
                        newSet = CTopologySet;
                        newSet.mMasses = massPros;
                        newSet.mRootMono = newUnit2;
                        newSet.mGapMono{b} = newUnit1; % Add gap unit to branch{b}.
                        newSet.mSources = obj.mCurrentBranches;
                        newSet.mReconstructor = obj;
                        newSet.mMinusH = flagMinusH(a);
                        
                        if obj.mFinalPeak
                            newSet.mType = 'T';
                        else
                            newSet.mType = 'B';
                        end
                        obj.insert_into_currentTopologySuperSet( newSet );
                    end
                end
                
                if ~obj.mFinalPeak && obj.mTryCIon % then try C-ion
                    theoreticalMass = obj.mDelta2.B2C.mass + [massCompensationLow; massCompensationHigh];
                    %theoreticalMass = [obj.mDelta2.B2C.mass * (1 - obj.mMassAccuracyPP) + massCompensationLow; ...
                    %                   obj.mDelta2.B2C.mass * (1 + obj.mMassAccuracyPP) + massCompensationHigh];
                    idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                                ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
                    flagMinusH = zeros(1, length(idx));
                    if obj.mCheckMinus2H
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                    end
                    if obj.mCheckMinusH
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))];
                    end
                    
                    for a = 1 : length(idx)
                        newUnit1 = obj.mDelta2.B2B.unit(idx(a), 1);
                        if newUnit1.mClassID == 2  % Rule: Fuc can only be a branch, not the end of a linear structure.
                            continue;
                        end
                        newUnit2 = obj.mDelta2.B2B.unit(idx(a), 2);
                        if obj.acceptable_monosaccharide( newUnit1.mClass ) == 0 || ...
                                obj.acceptable_monosaccharide( newUnit2.mClass ) == 0
                            continue;
                        end
                        
                        temp = flagMinusH(a)*CMass.H;
                        massPros = [obj.mCurrentTargetMass, ...
                                    max(theoreticalMass(1, idx(a)), obj.mCurrentTargetMassLow+temp), ...
                                    min(theoreticalMass(2, idx(a)), obj.mCurrentTargetMassHigh+temp)] - CMass.H2O; % Although it is C-ion, we use its B-ion mass.

                        newSet = CTopologySet;
                        newSet.mMasses = massPros;
                        newSet.mRootMono = newUnit2;
                        newSet.mMissingMono = newUnit1;
                        newSet.mSources = obj.mCurrentBranches;
                        newSet.mReconstructor = obj;
                        newSet.mMinusH = flagMinusH(a);
                        newSet.mType = 'C';
                        
                        obj.insert_into_currentTopologySuperSet( newSet );
                        % S: 2018/12/29. 2-branch-with-gap
                        result = 1;
                        % E: 2018/12/29. 2-branch-with-gap
                        
                        for b = 1 : 4 % push the mGapMono into one of branches
                            if isempty( obj.mCurrentBranches{b} )
                                break;
                            end
                            
                            newSet = CTopologySet;
                            newSet.mMasses = massPros;
                            newSet.mRootMono = newUnit2;
                            newSet.mGapMono{b} = newUnit1;
                            newSet.mSources = obj.mCurrentBranches;
                            newSet.mReconstructor = obj;
                            newSet.mMinusH = flagMinusH(a);
                            newSet.mType = 'C';
                            obj.insert_into_currentTopologySuperSet( newSet );
                        end
                    end
                end
            end
        end
        
        % S: 2018/12/29. 2-branch-with-gap
        function result = test_topologyset_2BranchWithGap( obj, massBound )
            linkageMassLoss = CMass.H2O + 3 * obj.mPermethylated * CMass.CH2;
            result = 0;
            massCompensationLow = massBound(1);
            massCompensationHigh = massBound(2);
            
            % Check if the branches together are already too heavy or still too light
            temp = obj.mCurrentTargetMass - massCompensationLow;
            if obj.mPermethylated
                if temp < 160
                    result = -1; return;
                end
                if obj.mCheckGap
                    if temp > 860 % max( obj.mDelta2.B2C )
                        return;
                    end
                elseif temp > 420 % max( obj.mDelta.B2C )
                    return;
                end
            else
                if temp < 105
                    result = -1; return;
                end
                if obj.mCheckGap
                    if temp > 660 % max( obj.mDelta2.B2C )
                        return;
                    end
                elseif temp > 335 % max( obj.mDelta.B2C )
                    return;
                end
            end
            
            % Allow 1 gap only when no-gap fails
            if obj.mCheckGap == 1
                if isempty( obj.mCurrentBranches{2} ) && ~isempty( obj.mCurrentBranches{1} ) % Rule: Fuc can only be a branch, in a linear substructure.
                    possible = 0;
                    for aTPS = obj.mCurrentBranches{1}.mTopologySets
                        if aTPS.mRootMono.mClassID ~= 2
                            possible = 1;
                            break;
                        end
                    end
                    if ~possible, return; end
                end
                
                % first try B-ions
                theoreticalMass = obj.mDelta2.B2B.mass + [massCompensationLow; massCompensationHigh];
                idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                            ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
                flagMinusH = zeros(1, length(idx));
                if ~obj.mFinalPeak
                    if obj.mCheckMinus2H
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                    end
                    if obj.mCheckMinusH
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))];
                    end
                end
                
                for a = 1 : length(idx)
                    newUnit1 = obj.mDelta2.B2B.unit(idx(a), 1);
                    if newUnit1.mClassID == 2  % Rule: Fuc can only be a branch, not in a linear substructure.
                        continue;
                    end
                    newUnit2 = obj.mDelta2.B2B.unit(idx(a), 2);

                    temp = flagMinusH(a)*CMass.H;
                    massPros = [obj.mCurrentTargetMass, ...
                                max(theoreticalMass(1, idx(a)), obj.mCurrentTargetMassLow + temp), ...
                                min(theoreticalMass(2, idx(a)), obj.mCurrentTargetMassHigh + temp)];

                    % There are 3 branches in obj.mCurrentBranches{b}. 
                    % So there are 3 possiblities: [1 2] 3, [1 3] 2, [2 3] 1
                    % [1 2] 3
                    gapBranch = CTopologySet;
                    gapBranch.mMasses = massPros - (newUnit2.mMass - linkageMassLoss) - ...
                                       (obj.mCurrentBranches{3}.mMassPros(2) - CMass.Proton);
                    gapBranch.mRootMono = newUnit1;
                    gapBranch.mSources = {obj.mCurrentBranches{1}, obj.mCurrentBranches{2}, [], []};
                    gapBranch.mReconstructor = obj;
                    gapBranch.mType = 'B';
                    gapSS = CTopologySuperSet;
                    gapSS.mType = gapBranch.mType;
                    gapSS.mMassPros = gapBranch.mMasses;
                    gapSS.mTopologySets = gapBranch;
                    gapSS.mReconstructor = obj;
                    
                    newSet = CTopologySet;
                    newSet.mMasses = massPros;
                    newSet.mRootMono = newUnit2;
                    newSet.mSources = {gapSS, obj.mCurrentBranches{3}, [], []};
                    newSet.mReconstructor = obj;
                    newSet.mMinusH = flagMinusH(a);
                    newSet.mTargetPeaks = obj.mCurrentPeakIdx;
                    
                    if obj.mFinalPeak
                        newSet.mType = 'T';
                    else
                        newSet.mType = 'B';
                    end
                    obj.insert_into_currentTopologySuperSet( newSet );
                    
                    % [1 3] 2
                    if obj.mCurrentBranches{2} ~= obj.mCurrentBranches{3}
                        gapBranch = CTopologySet;
                        gapBranch.mMasses = massPros - (newUnit2.mMass - linkageMassLoss) - ...
                                            (obj.mCurrentBranches{2}.mMassPros(2) - CMass.Proton);
                        gapBranch.mRootMono = newUnit1;
                        gapBranch.mSources = {obj.mCurrentBranches{1}, obj.mCurrentBranches{3}, [], []};
                        gapBranch.mReconstructor = obj;
                        gapBranch.mType = 'B';
                        gapSS = CTopologySuperSet;
                        gapSS.mType = gapBranch.mType;
                        gapSS.mMassPros = gapBranch.mMasses;
                        gapSS.mTopologySets = gapBranch;
                        gapSS.mReconstructor = obj;
                        
                        newSet = CTopologySet;
                        newSet.mMasses = massPros;
                        newSet.mRootMono = newUnit2;
                        newSet.mSources = {gapSS, obj.mCurrentBranches{2}, [], []};
                        newSet.mReconstructor = obj;
                        newSet.mMinusH = flagMinusH(a);
                        newSet.mTargetPeaks = obj.mCurrentPeakIdx;
                        
                        if obj.mFinalPeak
                            newSet.mType = 'T';
                        else
                            newSet.mType = 'B';
                        end
                        obj.insert_into_currentTopologySuperSet( newSet );
                        
                        if obj.mCurrentBranches{1} ~= obj.mCurrentBranches{2}
                            % [2 3] 1
                            gapBranch = CTopologySet;
                            gapBranch.mMasses = massPros - (newUnit2.mMass - linkageMassLoss) - ...
                                                (obj.mCurrentBranches{1}.mMassPros(2) - CMass.Proton);
                            gapBranch.mRootMono = newUnit1;
                            gapBranch.mSources = {obj.mCurrentBranches{2}, obj.mCurrentBranches{3}, [], []};
                            gapBranch.mReconstructor = obj;
                            gapBranch.mType = 'B';
                            gapSS = CTopologySuperSet;
                            gapSS.mType = gapBranch.mType;
                            gapSS.mMassPros = gapBranch.mMasses;
                            gapSS.mTopologySets = gapBranch;
                            gapSS.mReconstructor = obj;
                            
                            newSet = CTopologySet;
                            newSet.mMasses = massPros;
                            newSet.mRootMono = newUnit2;
                            newSet.mSources = {gapSS, obj.mCurrentBranches{1}, [], []};
                            newSet.mReconstructor = obj;
                            newSet.mMinusH = flagMinusH(a);
                            newSet.mTargetPeaks = obj.mCurrentPeakIdx;
                            
                            if obj.mFinalPeak
                                newSet.mType = 'T';
                            else
                                newSet.mType = 'B';
                            end
                            obj.insert_into_currentTopologySuperSet( newSet );
                        end
                    elseif obj.mCurrentBranches{1} ~= obj.mCurrentBranches{3}
                    % [2 3] 1
                        gapBranch = CTopologySet;
                        gapBranch.mMasses = massPros - (newUnit2.mMass - linkageMassLoss) - ...
                                            (obj.mCurrentBranches{1}.mMassPros(2) - CMass.Proton);
                        gapBranch.mRootMono = newUnit1;
                        gapBranch.mSources = {obj.mCurrentBranches{2}, obj.mCurrentBranches{3}, [], []};
                        gapBranch.mReconstructor = obj;
                        gapBranch.mType = 'B';
                        gapSS = CTopologySuperSet;
                        gapSS.mType = gapBranch.mType;
                        gapSS.mMassPros = gapBranch.mMasses;
                        gapSS.mTopologySets = gapBranch;
                        gapSS.mReconstructor = obj;
                        
                        newSet = CTopologySet;
                        newSet.mMasses = massPros;
                        newSet.mRootMono = newUnit2;
                        newSet.mSources = {gapSS, obj.mCurrentBranches{1}, [], []};
                        newSet.mReconstructor = obj;
                        newSet.mMinusH = flagMinusH(a);
                        newSet.mTargetPeaks = obj.mCurrentPeakIdx;
                        
                        if obj.mFinalPeak
                            newSet.mType = 'T';
                        else
                            newSet.mType = 'B';
                        end
                        obj.insert_into_currentTopologySuperSet( newSet );
                    end
                end
                
                if ~obj.mFinalPeak && obj.mTryCIon % then try C-ion
                    theoreticalMass = obj.mDelta2.B2C.mass + [massCompensationLow; massCompensationHigh];
                    idx = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh ) & ...
                                ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow ) );
                    flagMinusH = zeros(1, length(idx));
                    if obj.mCheckMinus2H
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H2 ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H2 ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))*2];
                    end
                    if obj.mCheckMinusH
                        temp = find( ( theoreticalMass(1,:) < obj.mCurrentTargetMassHigh + CMass.H ) & ...
                                     ( theoreticalMass(2,:) > obj.mCurrentTargetMassLow + CMass.H ) );
                        idx = [idx, temp];
                        flagMinusH = [flagMinusH, ones(1, length(temp))];
                    end
                    
                    for a = 1 : length(idx)
                        %Need to fill this part
                    end
                end
            end
        end
        % E: 2018/12/29. 2-branch-with-gap
        
        function insert_into_currentTopologySuperSet( obj, newSet )
        % Add a CTopologySet newSet to a CTopologySuperSet of the same mass and same type 
        % in obj.mTopologySuperSets. Create a new CTopologySuperSet if necessary.
            if isempty( newSet ), return; end
            
            if newSet.mType == 'B' || newSet.mType == 'T' % insert into B-ions
                if ( newSet.mMinusH == 0 )
                    if isempty( obj.mCurrentTopologySuperSetB )
                        obj.mCurrentTopologySuperSetB = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetB;
                        TSS.mType = newSet.mType;
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        if newSet.mType == 'B'
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B];
                        else
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T];
                        end
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetB;
                    end
                elseif ( newSet.mMinusH == 1 )
                    if isempty( obj.mCurrentTopologySuperSetB1 )
                        obj.mCurrentTopologySuperSetB1 = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetB1;
                        TSS.mType = newSet.mType;
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        if newSet.mType == 'B'
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B_H]; % 11 = B-ion '-H'
                        else
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T_H]; % 31 = T-ion '-H'
                        end
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetB1;
                    end
                elseif ( newSet.mMinusH == 2 )
                    if isempty( obj.mCurrentTopologySuperSetB2 )
                        obj.mCurrentTopologySuperSetB2 = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetB2;
                        TSS.mType = newSet.mType;
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        if newSet.mType == 'B'
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.B_2H]; % 12 = B-ion '-2H'
                        else
                            TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.T_2H]; % 31 = T-ion '-2H'
                        end
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetB2;
                    end
                end
            elseif newSet.mType == 'C' % insert into C-ions
                if ( newSet.mMinusH == 0 )
                    if isempty( obj.mCurrentTopologySuperSetC )
                        obj.mCurrentTopologySuperSetC = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetC;
                        TSS.mType = 'C';
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C];
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetC;
                    end
                elseif ( newSet.mMinusH == 1 )
                    if isempty( obj.mCurrentTopologySuperSetC1 )
                        obj.mCurrentTopologySuperSetC1 = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetC1;
                        TSS.mType = 'C';
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C_H];
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetC1;
                    end
                elseif ( newSet.mMinusH == 2 )
                    if isempty( obj.mCurrentTopologySuperSetC2 )
                        obj.mCurrentTopologySuperSetC2 = CTopologySuperSet;
                        TSS = obj.mCurrentTopologySuperSetC2;
                        TSS.mType = 'C';
                        TSS.mMassPros(1) = obj.mCurrentTargetMass;
                        TSS.mTargetPeaks = [obj.mCurrentPeakIdx; CPeakType.C_2H];
                        TSS.mReconstructor = obj;
                    else
                        TSS = obj.mCurrentTopologySuperSetC2;
                    end
                end
            end
            
            TSS.add_aTopologySet( newSet );            
        end

        function add_currentTSS_to_pool( obj ) 
        % add mCurrentTopologySuperSetB/B1/B2/C/C1/C2 to the search space for the peaks followed
            peak = obj.mPeaks(obj.mCurrentPeakIdx);
            if ~isempty( obj.mCurrentTopologySuperSetC2 )
                if isempty( obj.mTopologySuperSets )
                    obj.mTopologySuperSets = obj.mCurrentTopologySuperSetC2;
                    peak.mInferredSuperSet = obj.mCurrentTopologySuperSetC2;
                else
                    if obj.mCurrentTopologySuperSetC2.mMassPros(2) > obj.mTopologySuperSets(end).mMassPros(2)
                        obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetC2];
                        peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC2];
                    else
                        for insertIdx = length(obj.mTopologySuperSets) : -1 : 1
                            if obj.mTopologySuperSets(insertIdx).contains( obj.mCurrentTopologySuperSetC2 )
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mTopologySuperSets(insertIdx)];
                                obj.mTopologySuperSets(insertIdx).add_peak( obj.mCurrentPeakIdx, 22 );
                                break;
                            elseif obj.mCurrentTopologySuperSetC2.mMassPros(2) > obj.mTopologySuperSets(insertIdx).mMassPros(2)
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC2];
                                obj.mTopologySuperSets = [obj.mTopologySuperSets(1:insertIdx), obj.mCurrentTopologySuperSetC2, obj.mTopologySuperSets(insertIdx+1:end)];
                                break;
                            end
                        end
                    end
                end
            end
            if ~isempty( obj.mCurrentTopologySuperSetC1 )
                if isempty( obj.mTopologySuperSets )
                    obj.mTopologySuperSets = obj.mCurrentTopologySuperSetC1;
                    peak.mInferredSuperSet = obj.mCurrentTopologySuperSetC1;
                else
                    if obj.mCurrentTopologySuperSetC1.mMassPros(2) > obj.mTopologySuperSets(end).mMassPros(2)
                        obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetC1];
                        peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC1];
                    else
                        for insertIdx = length(obj.mTopologySuperSets) : -1 : 1
                            if obj.mTopologySuperSets(insertIdx).contains( obj.mCurrentTopologySuperSetC1 )
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mTopologySuperSets(insertIdx)];
                                obj.mTopologySuperSets(insertIdx).add_peak( obj.mCurrentPeakIdx, 12 );
                                break;
                            elseif obj.mCurrentTopologySuperSetC1.mMassPros(2) > obj.mTopologySuperSets(insertIdx).mMassPros(2)
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC1];
                                obj.mTopologySuperSets = [obj.mTopologySuperSets(1:insertIdx), obj.mCurrentTopologySuperSetC1, obj.mTopologySuperSets(insertIdx+1:end)];
                                break;
                            end
                        end
                    end
                end
            end
            if ~isempty( obj.mCurrentTopologySuperSetC )
                if isempty( obj.mTopologySuperSets )
                    obj.mTopologySuperSets = obj.mCurrentTopologySuperSetC;
                    peak.mInferredSuperSet = obj.mCurrentTopologySuperSetC;
                else
                    if obj.mCurrentTopologySuperSetC.mMassPros(2) > obj.mTopologySuperSets(end).mMassPros(2)
                        obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetC];
                        peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC];
                    else
                        for insertIdx = length(obj.mTopologySuperSets) : -1 : 1
                            if obj.mTopologySuperSets(insertIdx).contains( obj.mCurrentTopologySuperSetC )
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mTopologySuperSets(insertIdx)];
                                obj.mTopologySuperSets(insertIdx).add_peak( obj.mCurrentPeakIdx, 2 );
                                break;
                            elseif obj.mCurrentTopologySuperSetC.mMassPros(2) > obj.mTopologySuperSets(insertIdx).mMassPros(2)
                                obj.mTopologySuperSets = [obj.mTopologySuperSets(1:insertIdx), obj.mCurrentTopologySuperSetC, obj.mTopologySuperSets(insertIdx+1:end)];
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetC];
                                break;
                            end
                        end
                    end
                end
            end
            
            if ~isempty( obj.mCurrentTopologySuperSetB2 )
                if isempty( obj.mTopologySuperSets )
                    obj.mTopologySuperSets = obj.mCurrentTopologySuperSetB2;
                    peak.mInferredSuperSet = obj.mCurrentTopologySuperSetB2;
                else
                    if obj.mCurrentTopologySuperSetB2.mMassPros(2) > obj.mTopologySuperSets(end).mMassPros(2)
                        obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetB2];
                        peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetB2];
                    else
                        for insertIdx = length(obj.mTopologySuperSets) : -1 : 1
                            if obj.mTopologySuperSets(insertIdx).contains( obj.mCurrentTopologySuperSetB2 )
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mTopologySuperSets(insertIdx)];
                                obj.mTopologySuperSets(insertIdx).add_peak( obj.mCurrentPeakIdx, 21 );
                                break;
                            elseif obj.mCurrentTopologySuperSetB2.mMassPros(2) > obj.mTopologySuperSets(insertIdx).mMassPros(2)
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetB2];
                                obj.mTopologySuperSets = [obj.mTopologySuperSets(1:insertIdx), obj.mCurrentTopologySuperSetB2, obj.mTopologySuperSets(insertIdx+1:end)];
                                break;
                            end
                        end
                    end
                end
            end
            if ~isempty( obj.mCurrentTopologySuperSetB1 )
                if isempty( obj.mTopologySuperSets )
                    obj.mTopologySuperSets = obj.mCurrentTopologySuperSetB1;
                    peak.mInferredSuperSet = obj.mCurrentTopologySuperSetB1;
                else
                    if obj.mCurrentTopologySuperSetB1.mMassPros(2) > obj.mTopologySuperSets(end).mMassPros(2)
                        obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetB1];
                        peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetB1];
                    else
                        for insertIdx = length(obj.mTopologySuperSets) : -1 : 1
                            if obj.mTopologySuperSets(insertIdx).contains( obj.mCurrentTopologySuperSetB1 )
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mTopologySuperSets(insertIdx)];
                                obj.mTopologySuperSets(insertIdx).add_peak( obj.mCurrentPeakIdx, 11 );
                                break;
                            elseif obj.mCurrentTopologySuperSetB1.mMassPros(2) > obj.mTopologySuperSets(insertIdx).mMassPros(2)
                                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetB1];
                                obj.mTopologySuperSets = [obj.mTopologySuperSets(1:insertIdx), obj.mCurrentTopologySuperSetB1, obj.mTopologySuperSets(insertIdx+1:end)];
                                break;
                            end
                        end
                    end
                end
            end
            if ~isempty( obj.mCurrentTopologySuperSetB )
                obj.mTopologySuperSets = [obj.mTopologySuperSets, obj.mCurrentTopologySuperSetB];
                peak.mInferredSuperSet = [peak.mInferredSuperSet, obj.mCurrentTopologySuperSetB];
            end
        end
        
        function append_NLinkedRoot( obj )
            % target_mass_low = obj.mPeaks(end).mMass - obj.mMassAccuracyPPM;
            % target_mass_high = obj.mPeaks(end).mMass + obj.mMassAccuracyPPM;
            mass_error = obj.mMassAccuracyPP * obj.mPeaks(end).mMass;
            target_mass_low = obj.mPeaks(end).mMass - mass_error;
            target_mass_high = obj.mPeaks(end).mMass + mass_error;
            dMassBIonNoFuc = CMonosaccharideSet.get_mass( 'HexNAc', obj.mPermethylated ) * 2 - ...
                             (CMass.CH2 * 2 * obj.mPermethylated + CMass.H2O) * 2 + obj.mFinalPeakCompensation; %...
                             % obj.mReducingEndCompensation + CMass.CH2 * obj.mPermethylated + CMass.H2O;
            
            % check if there is a ion corresponding to Fuc
            superSet_Fuc = [];
            for k = 1 : length( obj.mTopologySuperSets )
                aSuperSet = obj.mTopologySuperSets(k);
                for cs = aSuperSet.mTopologySets
                    if cs.mRootMono.mClassID == 2 && isempty(cs.mSources{1}) && isempty(cs.mSources{1})
                        superSet_Fuc = aSuperSet;
                        break;
                    end
                end
                if ~isempty( superSet_Fuc )
                    break;
                end
            end
            if ~isempty( superSet_Fuc )
                dMassBIonWithFuc = dMassBIonNoFuc + CMonosaccharideSet.get_mass( 'Fuc', obj.mPermethylated ) - ...
                                   (CMass.CH2 * 2 * obj.mPermethylated + CMass.H2O);
            end
            
            num = length(obj.mTopologySuperSets);
            for cidx = num : -1 : 1
                aSuperSet = obj.mTopologySuperSets(cidx);
                tMassPros = aSuperSet.mMassPros + dMassBIonNoFuc;
                if tMassPros(3) < target_mass_low - 2 && isempty( superSet_Fuc ) % too light, no need to continue.
                    break;
                end
                
                if target_mass_high > tMassPros(2) && target_mass_low < tMassPros(3)
                    % create a 'aSuperSet - HexNAc - HexNAc' candidate set 
                    newTPSet = CTopologySet;
                    newTPSet.mType = 'T';
                    newTPSet.mMasses = tMassPros;
                    newTPSet.mMasses(1) = obj.mCurrentTargetMass; % Needed to enable insert_into_currentTopologySuperSet( newTPSet )
                    newTPSet.mRootMono = obj.mDelta.B2B.unit(5);
                    newTPSet.mMissingMono = obj.mDelta.B2B.unit(5);
                    newTPSet.mSources{1} = aSuperSet;
                    newTPSet.mReconstructor = obj;
                    obj.insert_into_currentTopologySuperSet( newTPSet )
                end
                
                if ~isempty( superSet_Fuc )
                    tMassPros = aSuperSet.mMassPros + dMassBIonWithFuc;
                    tMassPros(1) = obj.mCurrentTargetMass; % Needed to enable insert_into_currentTopologySuperSet( newTPSet )
                    
                    if target_mass_high > tMassPros(2) && target_mass_low < tMassPros(3)
                        % create a '[superSet_Fuc] [aSuperSet] HexNAc HexNAc' topology set
                        newTPSet = CTopologySet;
                        newTPSet.mType = 'T';
                        newTPSet.mMasses = tMassPros;
                        newTPSet.mRootMono = obj.mDelta.B2B.unit(5);
                        newTPSet.mMissingMono = obj.mDelta.B2B.unit(5);
                        newTPSet.mSources{1} = superSet_Fuc;
                        newTPSet.mSources{2} = aSuperSet;
                        newTPSet.mReconstructor = obj;
                        obj.insert_into_currentTopologySuperSet( newTPSet )

                        % create a '[superSet_Fuc] [aSuperSet HexNAc] HexNAc' topology set
                        newTPSet = CTopologySet;
                        newTPSet.mType = 'T';
                        newTPSet.mMasses = tMassPros;
                        newTPSet.mRootMono = obj.mDelta.B2B.unit(5);
                        newTPSet.mSources{1} = superSet_Fuc;
                        newTPSet.mSources{2} = aSuperSet;
                        newTPSet.mGapMono{2} = obj.mDelta.B2B.unit(5);
                        newTPSet.mReconstructor = obj;
                        obj.insert_into_currentTopologySuperSet( newTPSet )
                    end
                end
            end  % for cidx = 1 : length( obj.mTopologySuperSets )
        end
        
        function reconstruct_formulas( obj )
            % collect superSets needed to be reconstructed. This saves lots
            % of time for large glycans
            frontier = obj.mPeaks(end).mInferredSuperSet;
            buffer = frontier;
            while ~isempty( frontier )
                newFrontier = CTopologySuperSet.empty(0,0);
                for aSuperSet = frontier
                    for n = 1 : length( aSuperSet.mTopologySets )
                        aSet = aSuperSet.mTopologySets(n);
                        aSet.mTargetPeaks = aSuperSet.mTargetPeaks; % This is important because C-ion support is passed on to B-ion
                        
                        % clean -2H. Added 2020/06/21
                        tempLen = size(aSet.mTargetPeaks, 2);
                        flag = ones(1, tempLen);
                        for k = 1 : tempLen
                            if aSet.mTargetPeaks(2, k) == CPeakType.B_H || aSet.mTargetPeaks(2, k) == CPeakType.B_2H
                                if any( aSet.mTargetPeaks(2,:) == CPeakType.B )
                                    flag(k) = 0;
                                end
                            elseif aSet.mTargetPeaks(2, k) == CPeakType.C_H || aSet.mTargetPeaks(2, k) == CPeakType.C_2H
                                if any( aSet.mTargetPeaks(2,:) == CPeakType.C )
                                    flag(k) = 0;
                                end
                            end
                        end
                        aSet.mTargetPeaks = aSet.mTargetPeaks(:, flag>0);
                        
                        for k = 1 : 4
                            if ~isempty( aSet.mSources{k} ) && ~any( buffer == aSet.mSources{k} )
                                newFrontier = [newFrontier, aSet.mSources{k}];
                            else
                                break;
                            end
                        end
                    end
                end
                frontier = unique( newFrontier );
                buffer = [buffer, frontier];
            end
            buffer = unique( buffer );
            bufferMass = zeros(1, length(buffer));
            for k = 1 : length(buffer)
                bufferMass(k) = buffer(k).mMassPros(1);
            end
            [~, idx] = sort( bufferMass );
            buffer = buffer(idx);

%             buffer = obj.mTopologySuperSets;

            endFlank = ' ---';
            fprintf( '--- Superset ' );
            num = length(buffer);
            for k = 1 : num
                msg = [num2str(k), '/', num2str(num), endFlank];
                fprintf( '%s', msg );
                buffer(k).reconstruct_formulas();
                if k < num
                    for m = 1 : length( msg )
                        fprintf( '\b' );
                    end
                end
            end
            fprintf( '\nCGlycoDeNovo::reconstruct_formulas done!\n' );
            
            % dump illegal topologySets and topologySuperSets
            for peak = obj.mPeaks
                if isempty( peak.mInferredSuperSet )
                    continue;
                end
                
                peak.mInferredSuperSet = peak.mInferredSuperSet( [peak.mInferredSuperSet.mLegal] > 0 );
                if isempty( peak.mInferredSuperSet )
                    continue;
                end
                if length( peak.mInferredSuperSet ) > 1
                    formulas = [];
                    masses = [];
                    scores = [];
                    for superSet = peak.mInferredSuperSet
                        if ~isempty( superSet.mTopologies )
                            formulas = [formulas, {superSet.mTopologies.mFormula}];
                            masses = [masses, superSet.mTopologies.mMass];
                            for s = 1 : length(superSet.mTopologies)
                                scores(end+1) = length( superSet.mTopologies(s).mSupportPeaks );
                            end
                        end
                    end
                    [peak.mInferredFormulas, idx] = unique( formulas );
                    peak.mInferredMasses = masses(idx);
                    peak.mInferredScores = scores(idx);
                else %if ~isempty( peak.mInferredSuperSet )
                    superSet = peak.mInferredSuperSet;
                    if ~isempty( superSet.mTopologies )
                        peak.mInferredFormulas = {superSet.mTopologies.mFormula};
                        peak.mInferredMasses = [superSet.mTopologies.mMass];
                        peak.mInferredScores = zeros(1, length(superSet.mTopologies));
                        for s = 1 : length(superSet.mTopologies)
                            peak.mInferredScores(s) = length( superSet.mTopologies(s).mSupportPeaks );
                        end
                    end
                end
            end
            
            for aSuperSet = obj.mPeaks(end).mInferredSuperSet
                aSuperSet.sort_topologies_by_supports();
            end
        end % function reconstruct_formulas( obj )
        
        function rescore_reconstruction( obj )
            aPeak = obj.mPeaks(end);
            if ~isempty( aPeak.mInferredSuperSet )
                if isempty( obj.mPeakWeights ) || sum( obj.mPeakWeights < 1 ) == 0
                    return
                end
                
                formulas = {};
                masses = [];
                scores = [];
                for aSuperSet = aPeak.mInferredSuperSet
                    for aTP = aSuperSet.mTopologies
                        aTP.mScore = sum( obj.mPeakWeights( aTP.mPeaks ) );
                    end
                    aSuperSet.sort_topologies_by_supports();
                    formulas = [formulas, {aSuperSet.mTopologies.mFormula}];
                    masses = [masses, aSuperSet.mTopologies.mMass];
                    scores = [scores, aSuperSet.mTopologies.mScore];
                end
                [aPeak.mInferredFormulas, idx] = unique( formulas );
                aPeak.mInferredMasses = masses(idx);
                aPeak.mInferredScores = scores(idx);
                [~, idx] = sort( aPeak.mInferredScores, 'descend' );
                aPeak.mInferredScores = aPeak.mInferredScores(idx);
                aPeak.mInferredFormulas = aPeak.mInferredFormulas(idx);
                aPeak.mInferredMasses = aPeak.mInferredMasses(idx);
            end
        end
        
    end % methods
end