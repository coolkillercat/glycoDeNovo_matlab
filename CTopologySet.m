classdef CTopologySet < handle  % By Pengyu Hong @ Brandeis University
% A set of topologies that are of the same type, the same root and superset branches, and the "same" mass
    properties
        mType = '';
        mMasses = [0, 0, 0]; % [cannonical_mass, B-ion low_mass, B-ion high_mass];
        mTargetPeaks = -1;
        mLegal = 1;
        mMinusH = 0;
        
        mRootMono = CMonosaccharide.empty(0,0);
        mRootMonoClassID = 0;
        mMissingMono = CMonosaccharide.empty(0,0);
        mGapMono = {[], [], [], []};  % -> mGapMono -> mRootMono, one for each source
        mSources = {[], [], [], []};  % !!! The sources (CTopologySuperSet) are sorted from lightest to heaviest by this.mReconstructor

        mTopologyFormulas = {};
        mTopologies = CTopology.empty(0,0);
        mReconstructed = 0;  % 1 if it has been reconstructed from mono and sources
        
        mReconstructor = CGlycoDeNovo.empty(1,0);
    end
    
    properties (SetAccess = private)
        mCandidateNum = 0;
        mBranchTopologies = cell(1, 4);
    end
    
    methods
        function result = isequal( obj, topologySet )
            result = 0;
            
            if any( abs(obj.mMasses(2:3) - topologySet.mMasses(2:3) ) > 0.00001 )
                return;
            end
            
            if isempty( obj.mRootMono )
                if ~isempty( topologySet.mRootMono )
                    return;
                end
            elseif isempty( topologySet.mRootMono )
                return;
            elseif obj.mRootMono ~= topologySet.mRootMono
                return;
            end
            
            if isempty( obj.mGapMono )
                if ~isempty( topologySet.mGapMono )
                    return;
                end
            elseif isempty( topologySet.mGapMono )
                return;
            elseif ~isequal( obj.mGapMono, topologySet.mGapMono )
                return;
            end
            
            result = 1;
            for k = 1 : length( obj.mSources )
                if isempty( obj.mSources{k} )
                    if ~isempty( topologySet.mSources{k} )
                        result = 0;
                        return;
                    end
                elseif isempty( topologySet.mSources{k} )
                    result = 1;
                    return;
                elseif ~isequal( obj.mSources{k}, topologySet.mSources{k} )
                    result = 0;
                    return;
                end
            end
        end
        
        function reconstruct_formulas( obj )
            if obj.mReconstructed
                return;
            end
            for k = 1 : 4
                s = obj.mSources{k};
                if ~isempty(s) 
                    if s.mReconstructed == 0
                        s.reconstruct_formulas();
                    end
                else
                    break;
                end
            end
            
            obj.mTopologyFormulas = cell(1, 2000); % reserve some memory
            obj.mTopologies = CTopology.empty(0, 2000);
            obj.mCandidateNum = 0;
            obj.rec_formula( 1 );
            obj.mReconstructed = 1;
            obj.mLegal = (obj.mCandidateNum > 0);
            obj.mTopologyFormulas = obj.mTopologyFormulas(1:obj.mCandidateNum);
            obj.mTopologies = obj.mTopologies(1:obj.mCandidateNum);
        end
        
        function printtopology(obj)
            fprintf('%d',obj.mRootMonoClassID);
            for k = 1 : 4
                 s = obj.mSources{k};
                 if ~isempty(s)
                     s.printtopology();
                     fprintf('|');
                 end
            end
        end
    end
       
    methods (Access = private)
        function rec_formula( obj, sourceIdx )
            if (sourceIdx > 4) || isempty( obj.mSources{sourceIdx} )  % boundary met, create a CTopology
                compositionCountMerged = zeros(1, CMonosaccharideSet.cNumberMonosaccharideClasses);
                compositionCountMerged( obj.mRootMono.mClassID ) = 1;
                % S: 2018/12/29. 2-branch-with-gap
                if isempty(obj.mTargetPeaks)
                    peaks = [];
                else
                    peaks = obj.mTargetPeaks(1,:);
                end
                % E: 2018/12/29. 2-branch-with-gap
                
                % if obj.mReconstructor.mPermethylated
                if obj.mRootMono.mPermethylated
                    mass = obj.mRootMono.mMass - CMass.H2O - CMass.CH2 + CMass.Proton;
                else
                    mass = obj.mRootMono.mMass - CMass.H2O + CMass.Proton;
                end
                
                if obj.mType == 'T'
                    mass = mass + obj.mReconstructor.mFinalPeakCompensation;
                end
                
                if isempty( obj.mMissingMono )
                    formula = obj.mRootMono.mClass;
                else
                    formula = [obj.mMissingMono.mClass, ' ', obj.mRootMono.mClass];
                    compositionCountMerged( obj.mMissingMono.mClassID ) = compositionCountMerged( obj.mMissingMono.mClassID ) + 1;
                    if obj.mReconstructor.mPermethylated
                        mass = mass + obj.mMissingMono.mMass - CMass.H2O - CMass.CH2 - CMass.CH2;
                    else
                        mass = mass + obj.mMissingMono.mMass - CMass.H2O;
                    end
                end
                
                if sourceIdx > 1
                    % append branch formulas to the root
                    branchFormula = cell(1, sourceIdx-1);
                    for kk = 1 : sourceIdx-1
                        if isempty( obj.mMissingMono )
                            if ~CGlycoDeNovo.isLegalGlycosidicBond( obj.mBranchTopologies{kk}.mRootMonoClassID, obj.mRootMono.mClassID ) 
                                return;
                            end
                        else
                            if ~CGlycoDeNovo.isLegalGlycosidicBond( obj.mBranchTopologies{kk}.mRootMonoClassID, obj.mMissingMono.mClassID ) 
                                return;
                            end
                        end
                        
                        % if obj.mReconstructor.mPermethylated
                        if obj.mRootMono.mPermethylated
                            mass = mass + obj.mBranchTopologies{kk}.mMass - CMass.CH2 - CMass.Proton;
                        else
                            mass = mass + obj.mBranchTopologies{kk}.mMass - CMass.Proton;
                        end
                        
                        peaks = [peaks, obj.mBranchTopologies{kk}.mSupportPeaks];
                        compositionCountMerged = compositionCountMerged + obj.mBranchTopologies{kk}.mCompositionCount;
                        if isempty( obj.mGapMono{kk} )
                            branchFormula{kk} = obj.mBranchTopologies{kk}.mFormula;
                        else
                            if obj.mReconstructor.mPermethylated
                                mass = mass + obj.mGapMono{kk}.mMass - CMass.H2O - CMass.CH2 - CMass.CH2;
                            else
                                mass = mass + obj.mGapMono{kk}.mMass - CMass.H2O;
                            end
                            branchFormula{kk} = [obj.mBranchTopologies{kk}.mFormula, ' ', obj.mGapMono{kk}.mClass];
                            compositionCountMerged( obj.mGapMono{kk}.mClassID ) = compositionCountMerged( obj.mGapMono{kk}.mClassID ) + 1;
                        end
                    end

                    if ( obj.mMasses(2) - obj.mMasses(3) > -0.000001 ) % Want to test if obj.mMasses(2) == obj.mMasses(3). However, due to rounding error, sometimes obj.mMasses(2) may not be exactly equal to obj.mMasses(3)
                        if ~isempty(obj.mReconstructor) && ( abs( mass - obj.mMasses(2) ) > obj.mReconstructor.mMassAccuracyDalton )
                            return;
                        end
                    elseif (mass <= obj.mMasses(2)) || (mass >= obj.mMasses(3))
                        return;
                    end
                    
                    if sourceIdx == 2
                        formula = [branchFormula{1}, ' ', formula];
                    else
                        branchFormula = sort( branchFormula );
                        for kk = 1 : length( branchFormula )
                            formula = ['[', branchFormula{kk}, '] ', formula];
                        end
                    end
                end
                
                if isempty(obj.mReconstructor) || ~any( compositionCountMerged > obj.mReconstructor.mCompositionConstraint )
                    if ~isempty(obj.mReconstructor) && obj.mType == 'T' && obj.mReconstructor.mNLinked
                        if ~CGlycan.checkNLinked( formula )
                            return;
                        end
                    end
                    if any( strcmp( obj.mTopologyFormulas(1:obj.mCandidateNum), formula ) )
                        return;
                    end

                    newTP = CTopology;
                    newTP.mFormula = formula;
                    newTP.mRootMonoClassID = obj.mRootMono.mClassID;
                    newTP.mSupportPeaks = unique( peaks(1,:) );
                    newTP.mScore = length( newTP.mSupportPeaks );
                    newTP.mMass = mass;
                    newTP.mType = obj.mType;
                    newTP.mCompositionCount = compositionCountMerged;
                    
                    obj.mCandidateNum = obj.mCandidateNum + 1;
                    obj.mTopologies(obj.mCandidateNum) = newTP;
                    obj.mTopologyFormulas{ obj.mCandidateNum } = formula;
                end
            elseif ~isempty( obj.mSources{sourceIdx}.mTopologies ) % some turn out to be illegal during reconstruction
                for aTopology = obj.mSources{sourceIdx}.mTopologies
                    obj.mBranchTopologies{ sourceIdx } = aTopology; % iterate through all branch topologies
                    obj.rec_formula( sourceIdx + 1 );
                end
            end
        end
    end
end
