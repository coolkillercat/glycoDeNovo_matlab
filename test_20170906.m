datapath = 'D:\Projects\Glycomics\data\20170906\';
temp = dir( [datapath, '*.txt'] );
datafiles = cell(1, length(temp));
resultfiles = cell(1, length(temp));
for k = 1 : length(temp)
    datafiles{k} = [datapath, strrep(temp(k).name, '.txt', '')];
    resultfiles{k} = [datapath, 'kca.', strrep(temp(k).name, '.txt', '.mat')];
end

%%
fileIdx = 4;
kca = CKernelComponentAnalysis.process( datafiles{fileIdx} );

%%
save( resultfiles{fileIdx}, 'kca' );

%% batch
for fileIdx = 1 : 4
    kca = CKernelComponentAnalysis.process( datafiles{fileIdx} );
    save( resultfiles{fileIdx}, 'kca' );
end

%% Try several times to get the best estimates.
kca.load( datafiles{4} );
kca.analyze( 5 );
kca.compare_estimated_and_source( [], [150, 210]);

%%
figure; hold on;
plot(kca.mMassMap(:, 159)/max(kca.mMassMap(:, 159)) + kca.mMassMap(:, 142)/max(kca.mMassMap(:, 142)), 'k--', 'LineWidth', 2);
plot(kca.mMassMap(:, 142)/max(kca.mMassMap(:, 142)), 'r'); 
plot(kca.mMassMap(:, 159)/max(kca.mMassMap(:, 159)), 'b');
plot(kca.mMassMap(:, end)/max(kca.mMassMap(:, end)), 'Color', [0, 0.5, 0]); hold off;
legend( {'442.2049+472.2155', '442.2049', '472.2155', 'precursor' } );
xlim([310, 410]);

%% select a few peaks for presentation.
dp = kca.detect_diagnosis_peaks;

allpeaks = [];
for k = 1 : length(dp)
    allpeaks = [allpeaks, [dp{k}.peak]];
end

allpeaks = unique(allpeaks);
maxintensity = max( kca.mFilteredMassMap(:, allpeaks), 2 );
count = zeros(1, length(allpeaks));
for k = 1 : length(dp)
    temp = [dp{k}.peak];
    for t = temp
        idx = find( t == allpeaks );
        count(idx) = count(idx) + 1;
    end
end
sp = allpeaks(count == 1);
smi = maxintensity(count == 1);
temp = sort(smi, 'descend');
sp = sp( smi > temp(25));

%%
kca.compare_estimated_and_source(sp);
