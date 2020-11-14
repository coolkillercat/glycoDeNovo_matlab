classdef CIonProcessor < handle % By Pengyu Hong @ Brandeis University
    
    properties
        mIonClassifier = [];
        mMassAccuracy = 0.02;
    end
    
    methods (Static)
        function [zscores, intensities] = standardize_intensity( peaks, transform, peakIDs )
            if nargin < 2, transform = ''; end
            if nargin < 3, peakIDs = []; end
            
            intensities = [peaks.mIntensity];
            switch transform
                case 'log'
                    intensities = log(intensities);
                case 'sqrt'
                    intensities = sqrt(intensities);
            end
            
            % [m, std] = robust_mean_std( intensities( [peaks.mComplement] == 0 ) );
            if isempty( peakIDs )
                [m, std] = robust_mean_std( intensities(1:end-1) );
                zscores = (intensities - m)/std;
            else
                [m, std] = robust_mean_std( intensities(peakIDs) );
                zscores = (intensities - m)/std;
            end
        end
    end
    
    methods
        function obj = CIonProcessor( rebuildClassifier )
            if nargin < 1 || isempty(rebuildClassifier), rebuildClassifier = 0; end
            if isempty( obj.mIonClassifier )
                if ~rebuildClassifier && exist( 'ionclassifier.mat', 'file' )
                    load( 'ionclassifier.mat' );
                    obj.mIonClassifier = mIonClassifier;
                else
                    obj.mIonClassifier = CIonClassifier;
                end
            end
        end
        
        function features = convert2features( this, masses )
            features = this.mIonClassifier.convert2featuers(masses);
        end
        
        function setFeatures( this, massFeatures, massAccuracy )
            this.mMassFeatures = massFeatures;
            this.mMassAccuracy = massAccuracy;
            this.mNumFeatures = length(massFeatures);
        end
        
        function score_topologies( this, spectra, useIonCorrelator )
        % function score_topologies( this, spectra, useIonCorrelator )
        
            if nargin < 3, useIonCorrelator = 0; end
            
            downDelta = min( this.mIonClassifier.mMassFeatures );
            upDelta = max( this.mIonClassifier.mMassFeatures );

            for k = 1 : length(spectra)
                if isempty( spectra(k).mPeaks(end).mInferredSuperSet )
                    continue;
                end
    
                peaks = spectra(k).mPeaks;
                masses = [peaks.mMass];
                intensities = CIonProcessor.standardize_intensity( peaks, 'log' );
                TSS = spectra(k).mPeaks(end).mInferredSuperSet;
   
                for tp = TSS.mTopologies
                    weights = zeros(1, length(tp.mSupportPeaks));
                    
                    for m = 1 : length(tp.mSupportPeaks)
                        pidx = tp.mSupportPeaks(m);
                        peakMass = masses( pidx );
                        targetPeaks = peaks(pidx).mInferredSuperSet.mTargetPeaks;
                        
                        switch targetPeaks(2, targetPeaks(1,:) == pidx )
                            case 1
                                type = 'B';
                                flag = ( masses >= peakMass + downDelta ) & ( masses <= peakMass + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - peakMass; intensities(flag)];
                            case 2
                                type = 'C';
                                flag = ( masses >= peakMass + downDelta ) & ( masses <= peakMass + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - peakMass; intensities(flag)];
                            case 11
                                type = 'B';
                                flag = ( masses >= peakMass + CMass.H + downDelta ) & ( masses <= peakMass + CMass.H + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - (peakMass + CMass.H); intensities(flag)];
                            case 12
                                type = 'C';
                                flag = ( masses >= peakMass + CMass.H + downDelta ) & ( masses <= peakMass + CMass.H + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - (peakMass + CMass.H); intensities(flag)];
                            case 21
                                type = 'B';
                                flag = ( masses >= peakMass + CMass.H2 + downDelta ) & ( masses <= peakMass + CMass.H2 + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - (peakMass + CMass.H2); intensities(flag)];
                            case 22
                                type = 'C';
                                flag = ( masses >= peakMass + CMass.H2 + downDelta ) & ( masses <= peakMass + CMass.H2 + upDelta );
                                flag(pidx) = 0;
                                vicinity = [masses(flag) - (peakMass + CMass.H2); intensities(flag)];
                        end
                        vicinity = [[peakMass; intensities(pidx)], vicinity];
                        c = this.mIonClassifier.classify( vicinity, type, 2 );
                        weights(m) = 2*(c.overall - 0.5);
%                         switch type
%                             case 'B'
%                                 weights(m) = c.BY + c.BZ - 1;
%                             case 'C'
%                                 weights(m) = c.CY + c.CZ - 1;
%                         end
                    end
                    tp.mScore = sum(weights);
                end
                TSS.sort_topologies_by_supports();
                
                % break tied by correlation with intensities.
                if useIonCorrelator
                    tiedN = 0;
                    tp = TSS.mTopologies(1);
                    for m = 2 : length(TSS.mTopologies)
                        % if abs(tp.mScore - TSS.mTopologies(m).mScore) < 0.01
                        if TSS.mTopologies(m).mScore / tp.mScore > 0.8
                            tiedN = tiedN + 1;
                        else
                            break;
                        end
                    end
                    if tiedN > 1
                        zscores = CIonProcessor.standardize_intensity( peaks, 'sqrt' );
                        masses = [peaks.mMass];
                        
                        for m = 1 : tiedN
                            tp = TSS.mTopologies(m);
                            tpMass = masses( tp.mSupportPeaks(1:end-1) );
                            tpZ = zscores( tp.mSupportPeaks(1:end-1) );
                            
                            g = CGlycan( spectra(k).mO18, spectra(k).mPermethylated, spectra(k).mDeuterium, spectra(k).mReducedEnd );
                            g.parse( tp.mFormula );
                            
                            [Ys, Bs, ~, ~, Zs, Cs] = g.cleave( 'Proton' ); % [Ys, Bs, Xs, As, Zs, Cs] = cleave( obj, specU.mMetal );
                            simCount = zeros(1, length(tpMass));
                            
                            for ionM = [[Bs.mass], [Cs.mass], [Ys.mass], [Zs.mass]]
                                idxM = find( abs( tpMass - ionM ) < this.mMassAccuracy );
                                if ~isempty( idxM  )
                                    simCount(idxM) = simCount(idxM) + 1 / length(idxM);
                                end
                            end
                            tp.mScore = tp.mScore * max(corr( tpZ', simCount' ), 0.05);
                        end
                        
                        for m = tiedN+1 : length(TSS.mTopologies)
                            tp = TSS.mTopologies(m);
                            tp.mScore = tp.mScore * 0.05;
                        end
                        
                        TSS.sort_topologies_by_supports();
                    end
                end
            end
        end
        
        function setIonClassifierUseIntensity( this, useIntensity ) % useIntensity = 0 / 1
            this.mIonClassifier.mUseIntensity = useIntensity;
        end
        
        function trainData = train( this, spectra ) % B{k}.mass, B{k}.vicinity
            this.mIonClassifier.mMassAccuracy = this.mMassAccuracy;
            if nargout > 0
                trainData = this.mIonClassifier.train( spectra );
            else
                this.mIonClassifier.train( spectra );
            end
        end
    end
end