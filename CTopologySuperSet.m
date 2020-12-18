classdef CTopologySuperSet < handle  % By Pengyu Hong @ Brandeis University
% A set of topologies that are of the same type and the "same" mass
    properties
        mType = '';
        mMassPros = [0, Inf, 0]; % [peak_mass, low_mass, high_mass]
        mTargetPeaks = []; % mTargetPeaks(1,:) -- peak ID; mTargetPeaks(2,:) -- ion type
        mLegal = 1;
        
        mFormulas = {};
        mTopologies = CTopology.empty(0,0);
        mTopologySets = CTopologySet.empty(0,0);
        mLinkageScore;
        mReconstructed = 0;  % 1 if it has been reconstructed from mono and sources

        mReconstructor = CGlycoDeNovo.empty(1,0);
    end
    
    methods
        
        function add_aTopologySet( obj, newSet )
            if newSet.mType ~= obj.mType
                return;
            end
            
            % We don't need this?
            for k = 1 : length( obj.mTopologySets )
                oldSet = obj.mTopologySets(k);
                if oldSet.isequal( newSet )
                    return;
                end
            end
            
            obj.mTopologySets(end+1) = newSet;
            obj.mMassPros(2) = min( newSet.mMasses(2), obj.mMassPros(2) );
            obj.mMassPros(3) = max( newSet.mMasses(3), obj.mMassPros(3) );
            
            newSet.mTargetPeaks = obj.mTargetPeaks;
        end
        
        function add_peak( obj, peak, peakType )
            if ~any( obj.mTargetPeaks(1,:) == peak ) % Modified 2020/03/26
                obj.mTargetPeaks(:, end+1) = [peak; peakType];
            end
        end
        
        function result = contains( obj, tss )
            result = 0;
            
            if tss.mMassPros(2) < obj.mMassPros(2) - 0.0000001 || tss.mMassPros(3) > obj.mMassPros(3) + 0.0000001
                return;
            end

            for tssTS = tss.mTopologySets
                notfound = 1;
                for objTS = obj.mTopologySets
                    if objTS.isequal( tssTS )
                        notfound = 0;
                        break;
                    end
                end
                if notfound
                    return;
                end
            end
            
            result = 1;
        end
        
        function reconstruct_formulas( obj )
            if obj.mReconstructed
                return;
            end
            
            obj.mFormulas = {};
            obj.mTopologies = CTopology.empty(0,0);
            num = length( obj.mTopologySets );
            flag = zeros(1, num);
            for k = 1 : num
                tpSet = obj.mTopologySets(k);
                tpSet.reconstruct_formulas();
                if ~isempty( tpSet.mTopologies )
                    flag(k) = 1;
                    obj.mTopologies = [obj.mTopologies, tpSet.mTopologies];
                end
            end
            % remove illegal tpSets
            obj.mReconstructed = 1;
            obj.mTopologySets = obj.mTopologySets( flag > 0 );
            obj.mLegal = ~isempty( obj.mTopologies );
            
            % merge redundance
            formulas = {obj.mTopologies.mFormula};
            [obj.mFormulas, idx, backIdx] = unique( formulas );
            if length( formulas ) > length( obj.mFormulas )
                for k = 1 : length(idx)
                    allidx = find( backIdx == backIdx(idx(k)) );
                    if length( allidx ) > 1
                        aTP = obj.mTopologies( allidx(1) );
                        for m = 2 : length( allidx )
                            aTP.mSupportPeaks = [aTP.mSupportPeaks, obj.mTopologies(allidx(m)).mSupportPeaks];
                        end
                        aTP.mSupportPeaks = unique( aTP.mSupportPeaks );
                    end
                end
            end
            obj.mTopologies = obj.mTopologies(idx);
        end
        
function reconstruct_formulas2( obj )
            if obj.mReconstructed
                return;
            end
            
            obj.mFormulas = {};
            obj.mTopologies = CTopology.empty(0,0);
            num = length( obj.mTopologySets );
            flag = zeros(1, num);
            for k = 1 : num
                tpSet = obj.mTopologySets(k);
                tpSet.reconstruct_formulas();
                if ~isempty( tpSet.mTopologies )
                    flag(k) = 1;
                    obj.mTopologies = [obj.mTopologies, tpSet.mTopologies];
                end
            end
            % remove illegal tpSets
            obj.mReconstructed = 1;
            obj.mTopologySets = obj.mTopologySets( flag > 0 );
            obj.mLegal = ~isempty( obj.mTopologies );
end
        
function printtopology(obj)
    for k = 1 : length(obj.mTopologySets)
        obj.mTopologySets(k).printtopology();
    end
end
        
        function sort_topologies_by_score( obj )
            supports = [obj.mTopologies.mScore];
            [~, idx] = sort( supports, 'descend' );
            obj.mTopologies = obj.mTopologies(idx);
            obj.mFormulas = {obj.mTopologies.mFormula};
        end
        
        function sort_topologies_by_supports( obj )
            supports = zeros(1, length(obj.mTopologies));
            for k = 1 : length(obj.mTopologies)
                supports(k) = length(obj.mTopologies(k).mSupportPeaks);
            end
            [~, idx] = sort( supports, 'descend' );
            obj.mTopologies = obj.mTopologies(idx);
            obj.mFormulas = {obj.mTopologies.mFormula};
        end
    end
    
    methods (Access = private)
        function result = add_topologies( obj, topologies )
            if isempty( obj.mTopologies )
                obj.mTopologies = topologies;
                result = 1;
            else
                oldFormulas = { obj.mTopologies.formula };
                flag = zeros(1, length(topologies));
                for k = 1 : length(topologies)
                    idx = find( strcmp( topologies(k).formula, oldFormulas ) );
                    if isempty(idx)
                        flag(k) = 1;
                    else
                        obj.mTopologies(idx).peaks = unique( [obj.mTopologies(idx).peaks, topologies(k).peaks] );
                    end
                end
                obj.mTopologies = [obj.mTopologies, topologies(flag>0)];
                result = sum(flag) > 0;
            end
        end
    end    
end
