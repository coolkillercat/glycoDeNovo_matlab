% data sets
datasets = {'20170615.2', '20180822', '20181128', '20190116_FTMS', '20200304' };

if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/'];
else
    datapath = 'D:\Projects\Glycomics\data\';
end

filenames = {};
for k = 1 : length( datasets )
    adir = [datapath, datasets{k}, filesep];
    files = dir( [adir, '*.txt'] );
    for m = 1 : length(files)
        filenames{end+1} = [adir, files(m).name];
    end
end
resultdir = [datapath, 'annotation_202006', filesep];

%%
specSet = CSpectrum.empty(length(filenames), 0);
for f = 1 : length( filenames )
    disp( ['Loading ', filenames{f}] );
    specSet{f} = CSpectrum.load( filenames{f} );
    specSet{f}.protonate();
    specSet{f}.add_complementary_ions();
    specSet{f}.merge_peaks( 0.001 );
end
disp('Done loading raw data ...');


%% O-glycan
datasets = {'20200612' };

if ismac
    datapath = [getenv('HOME'), '/Documents/Projects/Glycomics/data/'];
else
    datapath = 'D:\Projects\Glycomics\data\';
end

filenames = {};
for k = 1 : length( datasets )
    adir = [datapath, datasets{k}, filesep];
    files = dir( [adir, '*.txt'] );
    for m = 1 : length(files)
        filenames{end+1} = [adir, files(m).name];
    end
end
resultdir = [datapath, 'annotation_202006', filesep];

%%
specSet = CSpectrum.empty(length(filenames), 0);
for f = 1 : length( filenames )
    disp( ['Loading ', filenames{f}] );
    specSet{f} = CSpectrum.load( filenames{f} );
    specSet{f}.protonate();
    specSet{f}.add_complementary_ions();
    specSet{f}.merge_peaks( 0.001 );
end
disp('Done loading raw data ...');

%% 
annotate_glycosic_bonds( specSet, resultdir );