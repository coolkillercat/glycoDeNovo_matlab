function toSpec = copy_reconstruction( sSpec, toSpec )

toMass = [toSpec.mPeaks.mMass];
map = zeros(1, sSpec.num);
massAccuPPM = 0.000005;
complementMassDelta = toSpec.mPrecursor * massAccuPPM;

for k = sSpec.num : -1 : 1 % Loop from the highest to the lowest to avoid mis-assignment,
    sPeak = sSpec.mPeaks(k);
    if isempty( sPeak.mInferredSuperSet )
        continue;
    end
    
    if k == sSpec.num && abs(sPeak.mMass - toSpec.mPeaks(end).mMass) < toSpec.mPrecursor * massAccuPPM
        midx = toSpec.num;
    else
        [mv, midx] = min( abs( toMass - sPeak.mMass ) );
        if mv > massAccuPPM * sPeak.mMass && mv > 0.0025
            if toSpec.mPeaks(midx).mIsComplement && mv > complementMassDelta
                disp( ['Missed ', num2str(k)] );
                continue;
            end
        end
    end
    
    map(k) = midx;
    
    % redirect. Important!
    flagNew = ones(1, length(sPeak.mInferredSuperSet));
    for m = 1 : length(sPeak.mInferredSuperSet)
        aTSS = sPeak.mInferredSuperSet(m);
        
        for toTSS = toSpec.mPeaks(midx).mInferredSuperSet
            if ~isempty(aTSS.mFormulas) && isempty(setdiff( aTSS.mFormulas, toTSS.mFormulas )) 
            % Discard it if offers nothing new.
            % Use "~isempty(aTSS.mFormulas)" to keep B/C ions without
            % formula due to fail reconstruction. These B/C ions may be
            % false positives.
                flagNew(m) = 0;
                break;
            end
        end
        
        if ~flagNew(m), continue; end
        
        idx = ( aTSS.mTargetPeaks(1,:) == k );
        if isempty(idx)
            assert(false);
        end
        aTSS.mTargetPeaks(1, idx) = midx;
    end
    
    % for debug purpose
    % for m = 1 : length( sPeak.mInferredSuperSet )
    %    aTSS = sPeak.mInferredSuperSet(m);
    %    for t = 1 : length(aTSS.mTopologies)
    %        tp = aTSS.mTopologies(t);
    %        if any( tp.mSupportPeaks > length(map) )
    %            disp( [k, m, t] );
    %        end
    %    end
    % end
    
    toSpec.mPeaks(midx).mInferredSuperSet = [toSpec.mPeaks(midx).mInferredSuperSet, sPeak.mInferredSuperSet(flagNew>0)];
    toSpec.mPeaks(midx).mInferredFormulas = unique( [toSpec.mPeaks(midx).mInferredFormulas, sPeak.mInferredFormulas] );
end

done =  CTopology.empty(0,0);
for k = 1 : sSpec.num
    if isempty(sSpec.mPeaks(k).mInferredSuperSet)
        continue;
    end
    
    for m = 1 : length(sSpec.mPeaks(k).mInferredSuperSet)
        aTSS = sSpec.mPeaks(k).mInferredSuperSet(m);
        for t = 1 : length(aTSS.mTopologies)
            tp = aTSS.mTopologies(t);
            if sum( done == tp )
                continue;
            end
            done(end+1) = tp;
            for p = 1 : length(tp.mSupportPeaks)
                try
                    tp.mSupportPeaks(p) = map( tp.mSupportPeaks(p) );
                catch
                    stop = 1;
                end
            end
        end
    end
end