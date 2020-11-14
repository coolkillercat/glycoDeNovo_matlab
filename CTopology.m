classdef CTopology < handle % By Pengyu Hong @ Brandeis University
    properties
        mType = '';
        mMass = 0;
        mMinusH = 0;

        mFormula = '';
        mSupportPeaks = [];
        mScore = 0;
        mCompositionCount = zeros(1, CMonosaccharideSet.cNumberMonosaccharideClasses );
        mRootMonoClassID = 0;
    end
    
    methods
        function result = isequal( obj, cand )
            if obj.mType == cand.mType && strcmp(obj.mFormula, cand.mFormula)
                result = 1;
            else
                result = 0;
            end
        end
    end
end