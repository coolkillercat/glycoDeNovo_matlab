datapath = 'D:\Projects\Glycomics\data\20180119\';
temp = dir( [datapath, '*.mzXML'] );
datafiles = cell(1, length(temp));
resultfiles = cell(1, length(temp));
for k = 1 : length(temp)
    datafiles{k} = [datapath, temp(k).name];
    resultfiles{k} = [datapath, 'kca.', strrep(temp(k).name, 'mzXML', 'mat')];
end

%%
fileIdx = 1;
kca = CKernelComponentAnalysis.process( datafiles{fileIdx} );
