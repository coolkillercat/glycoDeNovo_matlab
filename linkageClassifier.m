classdef linkageClassifier
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        mMassFeatures = [];
        mMassAccuracy = 0.005;
    end
    
    methods
        
    end
    
    methods(Static)
        
        function [dataVectors, dataSources] = prepare_training_data(traindata)
            dataVectors = containers.Map;
            dataSources = containers.Map;
            for k = 1 : length(traindata)
                peak = traindata(k);
                linkage = [peak.type, peak.RE, peak.NRE];
                if ~isKey(dataVectors, linkage)
                    dataVectors(linkage) = [peak.mLinkageVector];
                    dataSources(linkage) = [k];
                else
                    dataVectors(linkage) = [dataVectors(linkage); peak.mLinkageVector];
                    dataSources(linkage) = [dataSources(linkage), k];
                end
            end
        end
        
        function massFeatures = calculate_mass_feature(specSet, masses, massBound, massAccuracy)
            if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 30]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            h1 = waitbar(0, 'spectrum');
            h2 = waitbar(0, 'peak');
            massFeatures = merge_masses(masses', massAccuracy, 1);
            for s = 1 : length(specSet)
                spectra = specSet{s};
                waitbar(s/length(specSet), h1 , ['spectrum progress...', num2str(100*s/length(specSet)), '%']);
                for k = 1 : length(spectra.mPeaks)
                    if (mod(k + 1, 10) == 0) 
                        waitbar(k/length(spectra.mPeaks), h2, ['peak progress...', num2str(100*k/length(spectra.mPeaks)), '%'])
                    end
                    peak = spectra.mPeaks(k);
                    masslen = length(massFeatures);
                    peak.mLinkageVector = [peak.mRawMZ, zeros(1, 2 * masslen)];
                    featureAvailable = ones(1, length(massFeatures));
                    neighbour = spectra.mPeaks(k).mNeighbour;
                    for kk = 1 : length(neighbour(1,:))
                        [d, idx] = min( abs( massFeatures - neighbour(1,kk) ) );
                        if d < massAccuracy + 0.001 && featureAvailable(idx)
                            peak.mLinkageVector(idx + 1) = 1;
                            peak.mLinkageVector(idx + 1 + masslen) = spectra.mPeaks(neighbour(2,kk)).mZscore;
                            spectra.mPeaks(k).mFeature = [spectra.mPeaks(k).mFeature, massFeatures(idx)];
                            spectra.mPeaks(k).mFeatureIntensity = [spectra.mPeaks(k).mFeatureIntensity, spectra.mPeaks(neighbour(2,kk)).mZscore];
                        end
                    end
                end
            end
            delete(h1);
            delete(h2);
        end
        
        function masses = calculate_masses(spectra,orimasses, massBound, massAccuracy)
             if nargin < 3 || isempty(massBound), massBound = [CMass.H*1.5, 30]; end
            if nargin < 4 || isempty(massAccuracy), massAccuracy = 0.005; end
            context = {};
            peakMasses = [spectra.mPeaks.mRawMZ];
            neighbour = cell(length(peakMasses), 1);
            peakAvailable = ones(1, length(peakMasses));
            [uniqueMasses, ~, uID] = unique( floor( peakMasses * 1000 ) / 1000 );
            numUniqueMasses = length(uniqueMasses);
            for kk = 1 : numUniqueMasses
                  idx = abs(peakMasses - uniqueMasses(kk)) < massAccuracy;
                  uniqueMasses(kk) = mean(peakMasses(idx));
            end
            uniqueMassAvailable = ones(1, numUniqueMasses);
            for k = 1 : numUniqueMasses
            	ionM = uniqueMasses(k);
            	[minDiff, minIdx] = min( abs( peakMasses - ionM ) );
                	if minDiff < massAccuracy
                        peakAvailable( minIdx ) = 0;
                    	uniqueMassAvailable(k) = 0;
                    	temp = abs(peakMasses - ionM);
                    	idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), minIdx );
                    	context{end+1} = [peakMasses(idx)-ionM; idx];
               	end
            end
            for k = 1 : length(peakMasses)
                pMass = peakMasses(k);
                [minDiff, minIdx] = min( abs( pMass - uniqueMasses ) );
                temp = abs(peakMasses - pMass);
                idx = setdiff( find( (temp < massBound(2)) & (temp > massBound(1)) ), k );
                neighbour{k} = [peakMasses(idx) - pMass; idx];
                spectra.mPeaks(k).mNeighbour = neighbour{k};
                if minDiff < massAccuracy && uniqueMassAvailable(minIdx) % Use uniqueMassAvailable to avoid duplication
                	peakAvailable( minIdx ) = 0;
                	uniqueMassAvailable(minIdx) = 0;
                    context{end+1} = neighbour{k};
               end
            end
            masses = orimasses;
            for k = 1 : length(context)
                masses = [masses, context{k}(1,:)];
            end
            masses = round( masses * 1000 ) / 1000;
            masses = unique(masses);
        end
    end
    
end

