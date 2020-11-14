%% 
massBound = 10; 
massAccuracy = 0.005;
savePath = 'E:\Brandeis\data\annotation\small_old.txt';
%% This dataset use SNAP peak picking

clear filename;
dataset = {'annotation_20200612'};
files = [];
numFiles = 0;
for f = 1 : length(dataset)
    if ismac
        datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', dataset{f}, filesep];
    else
        datapath = ['E:\Brandeis\data\', dataset{f}, filesep];
    end
        files = dir( [datapath, '*.ann'] );
        numFiles = numFiles + length(files);
end
%%
specSet = CSpectrum.empty(length(numFiles), 0);
total = 0;
pass = 0;
for f = 1 : length(dataset)
    if ismac
        datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/', dataset{f}, filesep];
    else
        datapath = ['E:\Brandeis\data\', dataset{f}, filesep];
    end
        files = dir( [datapath, '*.ann'] );
        for ff = 1 : length( files )
            filename = [datapath, files(ff).name];
            disp( ['Loading ', filename] );
            specSet{ff} = CSpectrum.load_ann( filename );
            total = total + length(specSet{ff}.mPeaks);
        end
end
total = total * 2;
disp('Done loading raw data ...');

%% merge mass
h=waitbar(0,'merge masses progress');
for s = 1 : length(specSet)
    spec = specSet{s};
    old = [];
    masses = [spec.mPeaks.mRawMZ];
    context = cell(length(masses),1);
    for m = 1 : length(masses)
        pass = pass + 1;
        if (mod(pass, 100)==0)
            waitbar(pass/total, h , ['Calculating context spec ' num2str(s) ' peak ' num2str(m) ' ' num2str(100*pass/total) '%']);
        end
        idxes = find(abs(masses - masses(m)) < massBound);
        context{m} = [masses(idxes) - masses(m);idxes]; %
        old = [old, context{m}(1,:)];
    end
    massFeatures = merge_masses(old', massAccuracy, 1);
    for m = 1 : length(masses)
        pass = pass + 1;
        if (mod(pass, 100)==0)
            waitbar(pass/total, h , ['Extracting Feature spec ' num2str(s) ' peak ' num2str(m) ' ' num2str(100*pass/total) '%']);
        end
        for f = 1 : length(context{m}(1,:))
            [d, idx] = min( abs( massFeatures - context{m}(1,f) ) );
            if d < massAccuracy + 0.001
                spec.mPeaks(m).mFeature = [spec.mPeaks(m).mFeature, massFeatures(idx)];
                spec.mPeaks(m).mFeatureIntensity = [spec.mPeaks(m).mFeatureIntensity, spec.mPeaks(context{m}(2,f)).mZscore];
            end
        end
    end
end
delete(h);
%% output
%save(specSet)?
fid = fopen(savePath, 'w');
fprintf(fid, 'm/z\tz\tintensity\tpeaktype\tion-type\tRE\tNRE\tlinkage\tion-formula\tfeature\tfeatureintensity\n');
for s = 1 : length(specSet)
    spec = specSet{s};
    for p = 1:length(spec.mPeaks)
        peak = spec.mPeaks(p);
        com = '';
        if peak.mIsComplement
            com = 'com';
        end
        if strcmp(peak.type, 'O') == 0
            fprintf(fid, '%-f\t%-d\t%-f\t%-s\t%-s\t%-6s\t%-6s\t%-d\t%-15s\t%-s\t%-s\n', peak.mRawMZ, peak.mRawZ, peak.mZscore, com, peak.type, peak.RE, peak.NRE, peak.linkage, peak.mComment, vec2str(peak.mFeature), vec2str(peak.mFeatureIntensity));
        end
    end
end
fclose(fid);
