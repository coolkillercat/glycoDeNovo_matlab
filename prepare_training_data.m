function [trainVectors, trainData] = prepare_training_data( spectra, massAccuracy, useOriginalPeak, includeREM )

if nargin < 3, useOriginalPeak = 0; end
if nargin < 4, includeREM = 1; end

trainData.B = [];
trainData.C = [];
trainData.Y = [];
trainData.Z = [];
trainData.O = [];

% massLow = CMass.C + massAccuracy;
massLow = CMass.H * 1.5;
massHigh = 105;
glycanNames = {};
for s = 1 : length(spectra)
    specU = spectra(s);
    if strcmp(specU.mExperimentMethod, 'EED' ) == 0
        continue;
    end

    numbers = zeros(1, 4);
    peakMasses = [specU.mPeaks.mMass];
    peakComplements = [specU.mPeaks.mComplement] < 0;
    % peakIntensities = [peaks.mIntensity];
    flag = zeros(1, length(specU.mPeaks));
    for p = 1 : length(specU.mPeaks)-1
        if ~isempty(specU.mPeaks(p).mInferredFormulas)
            flag(p) = 1;
            
            if useOriginalPeak && specU.mPeaks(p).mComplement < 0
                flag(p) = 0;
                flag(-specU.mPeaks(p).mComplement) = 1;
            end
        end
    end
    peakIntensities = CIonProcessor.standardize_intensity( specU.mPeaks, 'log', find(flag == 1) );
    
    g = CGlycan( specU.mPermethylated, specU.mReducingEndModification );
    idx = strfind( specU.comment, '.' );
    name = strtrim( specU.comment(idx+1:end) );
    try
        g.parse( name );
    catch e
        disp( e );
        disp( s );
        continue;
    end
    glycanNames{end+1} = name;    
    [Ys, Bs, ~, ~, Zs, Cs] = g.cleave( 'Proton' ); % [Ys, Bs, Xs, As, Zs, Cs] = cleave( obj, specU.mMetal );
    allREM = {'', CMass.cReducingEndModification_O18, CMass.cReducingEndModification_Deuterium, ...
        CMass.cReducingEndModification_Reduced, CMass.cReducingEndModification_Aminopyridine, ...
        CMass.cReducingEndModification_PRAGS};
    REM = find( strcmp( specU.mReducingEndModification, allREM ) );
    
    bMasses = unique([Bs.mMass]); [~, idx, ~] = unique( round( bMasses * 1000 ) / 1000 ); bMasses = bMasses(unique(idx));
    cMasses = unique([Cs.mMass]); [~, idx, ~] = unique( round( cMasses * 1000 ) / 1000 ); cMasses = cMasses(unique(idx));
    yMasses = unique([Ys.mMass]); [~, idx, ~] = unique( round( yMasses * 1000 ) / 1000 ); yMasses = yMasses(unique(idx));
    zMasses = unique([Zs.mMass]); [~, idx, ~] = unique( round( zMasses * 1000 ) / 1000 ); zMasses = zMasses(unique(idx));
    available = ones(1, length(specU.mPeaks));
    
    availB = ones(1, length(bMasses));
    for k = 1 : length( bMasses )
        ionM = bMasses(k);
        [minDiff, idxM] = min( abs( peakMasses - ionM ) );
        if minDiff < massAccuracy
            available( idxM ) = 0;
            availB(k) = 0;
            temp = abs(peakMasses - ionM);
            idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), idxM );
            trainData.B(end+1).mass = peakMasses(idxM);
            trainData.B(end).peakID = idxM;
            trainData.B(end).complementPeakID = specU.mPeaks(idxM).mComplement;
            trainData.B(end).intensity = peakIntensities(idxM);
            trainData.B(end).vicinity = [peakMasses(idx)-ionM; peakIntensities(idx); peakComplements(idx)];
            trainData.B(end).glycanMass = g.mMass;
            trainData.B(end).glycanID = s;
            trainData.B(end).REM = REM;
            trainData.B(end).spectrum = spectra(s);
            numbers(1) = numbers(1) + 1;
        end
    end
    
    availC = ones(1, length(cMasses));
    for k = 1 : length(cMasses)
        ionM = cMasses(k);
        [minDiff, idxM] = min( abs( peakMasses - ionM ) );
        % if available( idxM ) && minDiff < massAccuracy
        if minDiff < massAccuracy
            available( idxM ) = 0;
            availC(k) = 0;
            temp = abs(peakMasses - ionM);
            idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), idxM );
            trainData.C(end+1).mass = peakMasses(idxM);
            trainData.C(end).peakID = idxM;
            trainData.C(end).complementPeakID = specU.mPeaks(idxM).mComplement;
            trainData.C(end).intensity = peakIntensities(idxM);
            trainData.C(end).vicinity = [peakMasses(idx)-ionM; peakIntensities(idx); peakComplements(idx)];
            trainData.C(end).glycanMass = g.mMass;
            trainData.C(end).glycanID = s;
            trainData.C(end).REM = REM;
            trainData.C(end).spectrum = spectra(s);
            numbers(2) = numbers(2) + 1;
        end
    end
    
    for k = 1 : length(specU.mPeaks)-1
        if (isempty( specU.mPeaks(k).mInferredFormulas ) || available( k ) == 0)
            continue; 
        end
        for TSS = specU.mPeaks(k).mInferredSuperSet
            type = TSS.mTargetPeaks(2, ( TSS.mTargetPeaks(1,:) == k ));
            switch type
                case {1, 2}
                    pMass = peakMasses(k);
                case {11, 21}
                    pMass = peakMasses(k) + CMass.H;
                case {12, 22}
                    pMass = peakMasses(k) + CMass.H2;
            end
            switch type
                case {1, 11, 12}
                    [minDiff, minIdx] = min( abs( pMass - bMasses ) );
                    if minDiff < massAccuracy && availB(minIdx)
                        available( minIdx ) = 0;
                        temp = abs(peakMasses - pMass);
                        idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), k );
                        trainData.B(end+1).mass = pMass;
                        trainData.B(end).peakID = k;
                        trainData.B(end).complementPeakID = specU.mPeaks(k).mComplement;
                        trainData.B(end).intensity = peakIntensities(k);
                        trainData.B(end).vicinity = [peakMasses(idx) - pMass; peakIntensities(idx); peakComplements(idx)];
                        trainData.B(end).glycanID = s;
                        trainData.B(end).glycanMass = g.mMass;
                        trainData.B(end).REM = REM;
                        trainData.B(end).spectrum = spectra(s);
                        % available(k) = 0;
                    end
                case {2, 21, 22}
                    [minDiff, minIdx] = min( abs( pMass - cMasses ) );
                    if minDiff < massAccuracy && availC(minIdx)
                        available( minIdx ) = 0;
                        temp = abs(peakMasses - pMass);
                        idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), k );
                        trainData.C(end+1).mass = pMass;
                        trainData.C(end).peakID = k;
                        trainData.C(end).complementPeakID = specU.mPeaks(k).mComplement;
                        trainData.C(end).intensity = peakIntensities(k);
                        trainData.C(end).vicinity = [peakMasses(idx) - pMass; peakIntensities(idx); peakComplements(idx)];
                        trainData.C(end).glycanID = s;
                        trainData.C(end).glycanMass = g.mMass;
                        trainData.C(end).REM = REM;
                        trainData.C(end).spectrum = spectra(s);
                        % available(k) = 0;
                    end
            end
        end
    end
        
    for ionM = yMasses
        [minDiff, idxM] = min( abs( peakMasses - ionM ) );
        % if available( idxM ) && minDiff < massAccuracy
        if minDiff < massAccuracy
            available( idxM ) = 0;
            temp = abs(peakMasses - ionM);
            idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), idxM );
            trainData.Y(end+1).mass = peakMasses(idxM);
            trainData.Y(end).peakID = idxM;
            trainData.Y(end).complementPeakID = specU.mPeaks(idxM).mComplement;
            trainData.Y(end).intensity = peakIntensities(idxM);
            trainData.Y(end).vicinity = [peakMasses(idx)-peakMasses(idxM); peakIntensities(idx); peakComplements(idx)];
            trainData.Y(end).glycanID = s;
            trainData.Y(end).glycanMass = g.mMass;
            trainData.Y(end).REM = REM;
            trainData.Y(end).spectrum = spectra(s);
            numbers(3) = numbers(3) + 1;
        end
    end
    if numbers(3) == 0
        disp( ['Missing Y: ', name]);
    end
    
    for ionM = zMasses
        [minDiff, idxM] = min( abs( peakMasses - ionM ) );
        % if available( idxM ) && minDiff < massAccuracy
        if minDiff < massAccuracy
            available( idxM ) = 0;
            temp = abs(peakMasses - ionM);
            idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), idxM );
            trainData.Z(end+1).mass = peakMasses(idxM);
            trainData.Z(end).peakID = idxM;
            trainData.Z(end).complementPeakID = specU.mPeaks(idxM).mComplement;
            trainData.Z(end).intensity = peakIntensities(idxM);
            trainData.Z(end).vicinity = [peakMasses(idx)-peakMasses(idxM); peakIntensities(idx); peakComplements(idx)];
            trainData.Z(end).glycanID = s;
            trainData.Z(end).glycanMass = g.mMass;
            trainData.Z(end).REM = REM;
            trainData.Z(end).spectrum = spectra(s);
            numbers(4) = numbers(4) + 1;
        end
    end
    if numbers(4) == 0
        disp( ['Missing Z: ', name]);
    end
    
%     for a = 1 : length(trainData.Z)-1
%         for b = a+1 : length(trainData.Z)
%             if trainData.Z{a}.peakID == trainData.Z{b}.peakID && trainData.Z{a}.glycanID == trainData.Z{b}.glycanID
%                 disp( ['Redundant ', num2str(s), ' ', num2str(a), ' ', num2str(b)] );
%             end
%         end
%     end

    disp( [name, ': ', sprintf( '%d ', numbers )] );
    
    for k = 1 : length(specU.mPeaks)-1
        if flag(k) && available( k )
            available(k) = 0;
            temp = abs( peakMasses - peakMasses(k) );
            idx = setdiff( find( (temp < massHigh) & (temp > massLow) ), k );
            trainData.O(end+1).mass = peakMasses(k);
            trainData.O(end).peakID = k;
            trainData.O(end).complementPeakID = specU.mPeaks(k).mComplement;
            trainData.O(end).intensity = peakIntensities(k);
            trainData.O(end).vicinity = [peakMasses(idx) - peakMasses(k); peakIntensities(idx); peakComplements(idx)];
            trainData.O(end).glycanID = s;
            trainData.O(end).glycanMass = g.mMass;
            trainData.O(end).REM = REM;
            trainData.O(end).spectrum = spectra(s);
        end
    end
end

% convert into vectors
trainVectors.B = [];
trainVectors.C = [];
trainVectors.Y = [];
trainVectors.Z = [];
trainVectors.O = [];

% build feature set
masses = [];
fields = { 'B', 'C', 'Y', 'Z', 'O' };
for f = 1 : length(fields)
    aField = fields{f};
    for k = 1 : length( trainData.(aField) )
        masses = [masses, trainData.(aField)(k).vicinity(1,:)];
    end
end
masses = round( masses * 1000 ) / 1000;
masses = unique(masses);
massFeatures = merge_masses(masses', massAccuracy, 1);
numFeatures = length( massFeatures );

for f = 1 : length( fields )
    aField = fields{f};
    num = length( trainData.(aField) );
    
    trainVectors.(aField) = zeros(num, numFeatures * 2 + 3);
    
    for k = 1 : num
        trainVectors.(aField)(k, end) = trainData.(aField)(k).REM;
        trainVectors.(aField)(k, end-1) = trainData.(aField)(k).glycanMass;
        trainVectors.(aField)(k, end-2) = trainData.(aField)(k).mass;
        for m = 1 : size( trainData.(aField)(k).vicinity, 2 )
            [d, idx] = min( abs( massFeatures - trainData.(aField)(k).vicinity(1, m) ) );
            if d < massAccuracy + 0.001
                trainVectors.(aField)(k, idx) = 1;
                trainVectors.(aField)(k, idx + numFeatures) = trainData.(aField)(k).vicinity(2, m);
            else
                disp( ['Missing mass features! ', aField, ': ', num2str(k)] );
            end
        end
    end
    
    if ~includeREM
        trainVectors.(aField) = trainVectors.(aField)(:, 1:end-1);
    end
end
trainVectors.massFeatures = massFeatures;
