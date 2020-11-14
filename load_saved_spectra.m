function spectra = load_saved_spectra( databatches )
if nargin < 1 || isempty( databatches )
    databatches = { '20120327', '20130530', '20130620', '20140904', '20141029', '20151121', '20160211', '20160410', '20170320', '20180426' };
end
if ~iscell( databatches )
    databatches = {databatches};
end

datapath = [strrep( pwd, 'matlab', 'data' ), filesep];

spectra = [];
dIdx = 1;
for b = 1 : length( databatches )
    files = dir( [datapath, databatches{b}, filesep, 'results', filesep, '*.mat'] );
    for f = 1 : length(files)
        datafilename = [datapath, databatches{b}, filesep, 'results', filesep, files(f).name];
        
        disp( datafilename );
        temp = load( datafilename );
        fieldnames = fields( temp );
        specU = temp.(fieldnames{1});
        specU.comment = strrep( specU.comment, char(9), ' ' );

        specU.filename = datafilename;
        spectra = [spectra, specU];
        dIdx = dIdx + 1;
    end
end
