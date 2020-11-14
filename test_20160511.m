%% Mixture of Glycans
amap = CSpectrumMap;
amap.load( '..\data\20160511\LNFP 1_3.23-55.txt' );
amap.source_decomposition();
amap.compare_data_with_source();

%%
specRange = amap.mSpectrumSet(1).mScanID : amap.mSpectrumSet(end).mScanID;
msheatmap( [minR:maxR]', sqrt(massmap(:, specRange)), 'SpecIdx', specRange, 'Markers', round(amap.mMassSet + CMass.Sodium - CMass.Proton) );


%% PCA decomposition
load( ['data', filesep, '20160511', filesep, 'unfiltered.mat'] );
[pcs, pcscores, recon, mu] = CSpectrumMap.mixture_decomposition_pca( unfiltered, 2 );

figure;
plot( [zeros(22,2); pcs] );
legend( {'Component 1', 'Component 2'} );
axis([23 55 -0.35 0.45]);
title( 'PCA Components' );

%%
figure;
for k = 1 : size( unfiltered, 2)
    subplot(2,4,k);
    plot( unfiltered(:,k), '--' ); hold on;
    plot( recon(:,k), 'k');
    plot( pcs(:,1) * pcscores(k,1) + mu'/2, 'r' );
    plot( pcs(:,2) * pcscoras(k,2) + mu'/2, 'b' );
    hold off;
    axis tight;
end
legend( {'Original', 'Estimation total', 'Estimation C1', 'Estimation C2'} );

%% Work direclty on the mzXML file (not successful so far)
mzxml_struct = mzxmlread( ['..', filesep, 'data', filesep, '20160511', filesep, 'lnfp1_3 cid.mzXML'] );
[peaks, ~] = mzxml2peaks(mzxml_struct, 'Levels', mzxml_struct.scan(1).msLevel);

%% Remove those heavier than the precursor
precursor = mzxml_struct.scan(1).precursorMz.value;
massAccuracy = 0.2;
for s = 1 : length( peaks )
    flag = (peaks{s}(:,1) < precursor + massAccuracy) & (peaks{s}(:,1) > 100);
    peaks{s} = peaks{s}(flag, :);
end

%% Correct baseline
for s = 1 : length( peaks )
    temp = msbackadj( peaks{s}(:,1), peaks{s}(:,2));
    temp(temp < 0)= 0;
    peaks{s}(:,2) = temp;
end

%% Remove trivial spectra
precursorIntensity = zeros(1, length(peaks));
for s = 1 : length(peaks)
    sidx = max( find( peaks{s}(:,1) < precursor - massAccuracy*2 ) );
    precursorIntensity(s) = max(peaks{s}(sidx:end,2));
end
precursorIntensity = log2(precursorIntensity);
intensityThreshold = median( precursorIntensity );

[~, midx] = max( precursorIntensity );
sidx = midx-1;
while sidx > 1 && (precursorIntensity(sidx) > intensityThreshold)
    sidx = sidx - 1;
end
sidx = max( [0, sidx-3] );

eidx = midx + 1;
while (eidx < length(peaks)) && (precursorIntensity(eidx) > intensityThreshold)
    eidx = eidx + 1;
end
eidx = min( [length(peaks), eidx+3] );

peaks = peaks(sidx:eidx);

%% Filter out background noise
filteredPeaks = peaks;
mv = zeros(1, length(filteredPeaks));
sv = zeros(1, length(filteredPeaks));
for s = 1 : length(filteredPeaks)
    mv(s) = mean(filteredPeaks{s}(:,2));
    sv(s) = std(filteredPeaks{s}(:,2));
end
mv = median(mv);
sv = median(sv);
for s = 1 : length(filteredPeaks)
    filteredPeaks{s}(filteredPeaks{s}(:,2) < mv+2*sv, 2) = 0;
end

%% Create MassMap
minMass = Inf; maxMass = -Inf;
for s = 1 : length(filteredPeaks)
    minMass = min( minMass, min(filteredPeaks{s}(:,1)) );
    maxMass = max( maxMass, max(filteredPeaks{s}(:,1)) );
end
masses = minMass : 0.1 : maxMass;

massMap = zeros(length(filteredPeaks), length(masses));
for s = 1 : length(filteredPeaks)
    intensity = zeros(1, length(masses));
    for n = 1 : size(filteredPeaks{s},1)
        temp = abs(filteredPeaks{s}(n,1) - masses);
        idxes = find(temp <= 1);
        weights = exp(-temp(idxes));
        intensity(idxes) = intensity(idxes) + filteredPeaks{s}(n,2) * (weights / sum(weights)) ;
    end
    massMap(s,:) = intensity;
end
