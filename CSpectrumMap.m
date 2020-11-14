classdef CSpectrumMap < handle
    properties
        mPrecursor = -1 % mz of the protonated precursor
        mSpectrumSet = CSpectrum.empty(0,0);
        
        mDeuterium = 0;
        mExperimentMethod = '';
        mMetal = '';
        mNLinked = 0;
        mO18 = 0;
        mPermethylated = 0;
        mProtonated = 0;
        mReducedEnd = 0;
        
        mMassMap = [];
        mMassSet = [];
        mSources = [];
        mSourceWeights = [];
        mEstimation = [];
        
        mComment = '';
    end
    
    methods
        function clear( obj )
            obj.mPrecursor = -1;
            obj.mSpectrumSet = CSpectrum.empty(0,0);
            obj.mDeuterium = 0;
            obj.mExperimentMethod = '';
            obj.mMetal = '';
            obj.mNLinked = 0;
            obj.mO18 = 0;
            obj.mPermethylated = 0;
            obj.mProtonated = 0;
            obj.mReducedEnd = 0;
            obj.mComment = '';
            obj.mMassMap = [];
            obj.mMassSet = [];
            obj.mSources = [];
            obj.mSourceWeights = [];
            obj.mEstimation = [];
        end
        
        function result = copy( obj )
            num = length(obj);
            result = CSpectrumMap.empty(0, num);
            if num == 0, return; end
            
            fnames = fields( obj(1) );
            for k = 1 : num
                result(k) = CSpectrumMap;
                for m = 1 : length(fnames)
                    if strcmp( fnames{m}, 'mSpectrumSet' )
                        result(k).mSpectrumSet = obj(k).mSpectrumSet.copy;
                    else
                        result(k).(fnames{m}) = obj(k).(fnames{m});
                    end
                end
            end
        end
        
        function [sources, approximated] = source_decomposition(obj, numSources)
            if nargin < 2 || isempty(numSources), numSources = 2; end
            
            [obj.mMassMap, obj.mMassSet] = obj.generate_compressed_map();
            [params, model] = CSpectrumMap.mixture_decomposition_gaussian_kernel( sqrt(obj.mMassMap), numSources );
            sources = [];
            for k = 1 : numSources
                sources.mean(k) = params( k*numSources-1 );
                sources.std(k) = params( k*numSources );
            end
            sources.model = model;
            obj.mSources = sources;
            obj.mSourceWeights = reshape( params(5:end), numSources, size(obj.mMassMap,2) )';
            [~, approximated] = model(params);
            obj.mEstimation = approximated .^ 2;
        end
        
        function compare_data_with_source( obj, idxes )
            if isempty( obj.mEstimation ) || isempty( obj.mMassMap )
                return;
            end
            
            if nargin < 2 || isempty( idxes )
                idxes = 1 : size( obj.mMassMap, 2 );
            end
            
            colors = 'rbmg';
            ori = obj.mMassMap(:, idxes);
            est = obj.mEstimation(:, idxes);
            selectedMap = obj.mMassSet(idxes);
            figure;
            x = 1 : size(ori, 1);
            num = size( ori, 2 );
            temp = sqrt(num / 12);
            rows = floor( temp * 3 );
            cols = ceil( num/rows );
            
            if obj.mSpectrumSet(1).mScanID > 0
                scanIdx = zeros(1, length(obj.mSpectrumSet));
                for k = 1 : length(scanIdx)
                    scanIdx(k) = obj.mSpectrumSet(k).mScanID;
                end
            else
                scanIdx = 1 : length(obj.mSpectrumSet);
            end

            legends= { 'Original', 'Est. Total' };
            for c = 1 : length(obj.mSources.mean)
                legends{end+1} = ['Est. Component ', num2str(c)];
            end
            for k = 1 : num
                subplot(rows, cols, k);
                plot( scanIdx, ori(:, k), '--', 'LineWidth', 1  ); hold on;
                plot( scanIdx, est(:, k), 'k' );
                cidx = 1;
                for c = 1 : length(obj.mSources.mean)
                    temp = obj.mSourceWeights(k,c) / obj.mSources.std(c) * exp(- (x-obj.mSources.mean(c)) .^ 2 / (2 * obj.mSources.std(c)^2) );
                    temp = temp .^ 2;
                    plot( scanIdx, temp, colors(cidx) );
                    cidx = cidx + 1;
                    if cidx > length(colors), cidx = length(colors); end
                end
                hold off;
                axis tight;
                title( ['m/z ', num2str(selectedMap(k), '%.3f')] );
                if k == 1
                    legend( legends, 'FontSize', 9 );
                end
            end
        end
        
        function [massMap, allMasses] = generate_compressed_map(obj)
            if obj.mProtonated
                allMasses = [];
                for aSpectrum = obj.mSpectrumSet
                    allMasses = [allMasses, [aSpectrum.mPeaks.mMass]];
                end
                allMasses = unique( allMasses )';
                
                massMap = zeros(length(obj.mSpectrumSet), length(allMasses));
                for k = 1 : length(obj.mSpectrumSet)
                    aSpectrum = obj.mSpectrumSet(k);
                    massSet = [aSpectrum.mPeaks.mMass];
                    for m = 1 : length(massSet)
                        idx = find( abs(massSet(m) - allMasses) < 0.0001 );
                        if isempty(idx)
                            throw( MException( 'CSpectrumMap::generate_compressed_map', 'Can not find mass ', double2str(massSet(m)) ) );
                        end
                        massMap(k, idx) = aSpectrum.mPeaks(m).mIntensity;
                    end
                end
            else
                allMasses = [];
                for aSpectrum = obj.mSpectrumSet
                    allMasses = [allMasses, [aSpectrum.mPeaks.mRawMZ]];
                end
                allMasses = unique( allMasses )';
                
                massMap = zeros(length(obj.mSpectrumSet), length(allMasses));
                for k = 1 : length(obj.mSpectrumSet)
                    aSpectrum = obj.mSpectrumSet(k);
                    massSet = [aSpectrum.mPeaks.mRawMZ];
                    for m = 1 : length(massSet)
                        idx = find( abs(massSet(m) - allMasses) < 0.0001 );
                        if isempty(idx)
                            throw( MException( 'CSpectrumMap::generate_compressed_map', 'Can not find mass ', double2str(massSet(m)) ) );
                        end
                        massMap(k, idx) = aSpectrum.mPeaks(m).mIntensity;
                    end
                end
            end
        end
        
        function merge_peaks( obj, threshold )
            if nargin < 2, threshold = 0.02; end
            
            allMasses = [];
            for aSpectrum = obj.mSpectrumSet
                allMasses = [allMasses, [aSpectrum.mPeaks.mMass]];
            end
            allMasses = unique( allMasses )';
            
            tree = linkage( allMasses, 'single' );
            c = cluster( tree, 'cutoff', threshold, 'criterion', 'distance' );
            newLen = max(c);
            newMasses = zeros(newLen, 1);
            
            for k = 1 : newLen
                tempM = allMasses( c == k );
                newMasses(k) = mean( tempM );
            end
           
            for k = 1 : length( obj.mSpectrumSet )
                flag = zeros(1, length(newMasses));
                intensity = zeros(1, length(newMasses));
                for m = 1 : obj.mSpectrumSet(k).num
                    [~, idx] = min( abs(newMasses - obj.mSpectrumSet(k).mPeaks(m).mMass ) );
                    flag(idx) = 1;
                    intensity(idx) = intensity(idx) + obj.mSpectrumSet(k).mPeaks(m).mIntensity;
                end
                flag = flag > 0;
                obj.mSpectrumSet(k) = CSpectrum.create_spectrum( newMasses(flag), intensity(flag) );
            end
        end
        
        function protonate( obj, considerOtherCharge, massAccuracy )
            if nargin < 2, considerOtherCharge = []; end
            if nargin < 3, massAccuracy = 0.01; end
            for aSpec = obj.mSpectrumSet
                aSpec.protonate( considerOtherCharge, massAccuracy );
            end
            obj.mProtonated = 1;
        end
        
        function load( obj, filename )
            obj.clear();
            
            fid = fopen( filename, 'r' );
            if fid == -1, disp( ['no such file: ', filename] ); return; end
            
            adductMetal = '';
            experimentMethod = '';
            O18Label = 0;
            permethylation = 0;
            reduced = 0;
            nlinked = 0;
            precursorMZ = -1;
            deu = 0;
            scanIDs = [];
            comment = '';
            
            aline = fgetl( fid );
            firstLine = 1;
            while str_startswith( aline, '#' ) || isempty(strtrim(aline))
                if str_startswith( aline, '# Metal:' )
                    adductMetal = strtrim( aline(9:end) );    
                    if strcmp( adductMetal, 'H' )
                        adductMetal = 'Proton';
                    end
                elseif str_startswith( aline, '# Method:' )
                    experimentMethod = strtrim( aline(10:end) );
                elseif str_startswith( aline, '# O18' )
                    O18Label = 1;
                elseif str_startswith( aline, '# Permethylated' )
                    permethylation = 1;
                elseif str_startswith( aline, '# Reduced' )
                    reduced = 1;
                elseif str_startswith( aline, '# NLinked' )
                    nlinked = 1;
                elseif str_startswith( aline, '# Deuterium' )
                    deu = 1;
                elseif str_startswith( aline, '# Precursor:' )
                    temp = strtrim( aline(13:end) );
                    fields = strsplit( temp, ';' );
                    for k = 1 : length(fields)
                        fields{k} = strtrim( fields{k} );
                        if str_startswith( fields{k}, 'H+' )
                            precursorMZ = str2double( fields{k}(3:end) );
                            break;
                        end
                    end
                elseif str_startswith( aline, '# Scans' )
                    aline = aline(9:end);
                    fields = strsplit( aline, '-' );
                    scanStart = str2double( fields{1} );
                    scanEnd = str2double( fields{2} );
                elseif str_startswith( aline, '# M/Z' )
                    fields = strsplit( aline, char(9) );
                    fields = fields(2:end);
                    formatSpec = '%f';
                    if length( fields ) == (scanEnd - scanStart) + 1
                        scanIDs = zeros(1, length(fields));
                        for k = 1 : length(scanIDs)
                            scanIDs(k) = str2double( fields{k} );
                            formatSpec = [formatSpec, '%f'];
                        end
                    else
                        scanIDs = zeros(1, length(fields)/2);
                        for k = 1 : length(scanIDs)
                            scanIDs(k) = str2double( fields{2*k-1} );
                            formatSpec = [formatSpec, '%f%f'];
                        end
                    end
                    break;
                elseif firstLine
                    comment = strtrim( aline(2:end) );
                    firstLine = 0;
                end
                
                aline = fgetl( fid );
            end
            
            data = textscan( fid, formatSpec );
            fclose( fid );
            
            if precursorMZ <= 0
                precursorMZ = max( data{1} );
            end
            
            for k = 1 : length(scanIDs)
                if length(scanIDs) == length(data)-1
                    aSpectrum = CSpectrum.create_spectrum( data{1}, data{k+1}, ones(length(data{1}),1) );
                else
                    aSpectrum = CSpectrum.create_spectrum( data{1}, data{2*k}, data{2*k+1} );
                end
                aSpectrum.mScanID = scanIDs(k);
                aSpectrum.mComment = comment;
                aSpectrum.mMetal = adductMetal;
                aSpectrum.mO18 = O18Label;
                aSpectrum.mReducedEnd = reduced;
                aSpectrum.mPermethylated = permethylation;
                aSpectrum.mExperimentMethod = experimentMethod;
                aSpectrum.mNLinked = nlinked;
                aSpectrum.mDeuterium = deu;
                aSpectrum.mPrecursor = precursorMZ;
                
                obj.mSpectrumSet(k) = aSpectrum;
            end
            
            obj.mComment = comment;
            obj.mMetal = adductMetal;
            obj.mO18 = O18Label;
            obj.mReducedEnd = reduced;
            obj.mPermethylated = permethylation;
            obj.mExperimentMethod = experimentMethod;
            obj.mNLinked = nlinked;
            obj.mDeuterium = deu;
            obj.mPrecursor = precursorMZ;
            obj.mComment = comment;
        end
        
        function load_xls( obj, filename )
            obj.clear();
            
            [~, sheets] = xlsfinfo(filename);            
            for k = 1 : length( sheets )
                if strcmpi( sheets{k}, 'annotation' )
                    [~, txt] = xlsread( filename, sheets{k} );
                    for m = 1 : length(txt)
                        aline = strtrim( txt{m} );
                        if str_startswith( aline, 'Metal:' )
                            chargeMetal = strtrim( aline(7:end) );
                            if strcmp( chargeMetal, 'H' ) || strcmp( chargeMetal, 'Proton' )
                                obj.mMetal = 'Proton';
                                obj.mProtonated = 1;
                            else
                                obj.mMetal = chargeMetal;
                            end
                        elseif str_startswith( aline, 'Method:' )
                            obj.mExperimentMethod = strtrim( txt{m}(8:end) );
                        elseif str_startswith( aline, 'O18' )
                            obj.mO18 = 1;
                        elseif str_startswith( aline, 'Permethylated' )
                            obj.mPermethylated = 1;
                        elseif str_startswith( aline, 'Reduced' )
                            obj.mReducedEnd = 1;
                        elseif str_startswith( aline, 'NLinked' )
                            obj.mNLinked = 1;
                        elseif str_startswith( aline, 'Deuterium' )
                            obj.mDeuterium = 1;
                        elseif str_startswith( aline, 'Precursor:' )
                            temp = strtrim( aline(11:end) );
                            fields = strsplit( temp, ';' );
                            for f = 1 : length(fields)
                                fields{f} = strtrim( fields{f} );
                                if str_startswith( fields{f}, 'H+' )
                                    obj.mPrecursor = str2double( fields{k}(3:end) );
                                    break;
                                end
                            end
                        elseif m == 1
                            obj.mComment = strtrim( aline(2:end) );
                        end
                    end
                else
                    data = xlsread( filename, sheets{k} );
                    aSpectrum = CSpectrum;
                    aSpectrum.mPrecursor = obj.mPrecursor; % mz of the protonated precursor
                    aSpectrum.mDeuterium = obj.mDeuterium;
                    aSpectrum.mExperimentMethod = obj.mExperimentMethod;
                    aSpectrum.mMetal = obj.mMetal;
                    aSpectrum.mNLinked = obj.mNLinked;
                    aSpectrum.mO18 = obj.mO18;
                    aSpectrum.mPermethylated = obj.mPermethylated;
                    aSpectrum.mProtonated = obj.mProtonated;
                    aSpectrum.mReducedEnd = obj.mReducedEnd;

                    aSpectrum.mPeaks = CPeak.empty(size(data, 1), 0);
                    for d = 1 : size(data, 1)
                        aSpectrum.mPeaks(d) = CPeak();
                        aSpectrum.mPeaks(d).mSpectrum = aSpectrum;
                        aSpectrum.mPeaks(d).mRawMZ = data(d,1);
                        aSpectrum.mPeaks(d).mRawZ = data(d,2);
                        aSpectrum.mPeaks(d).mIntensity = data(d,3);
                        if aSpectrum.mProtonated
                            aSpectrum.mPeaks(k).mMass = data(d,1);
                        end
                    end
                    
                    % in case the precursor is not included in the peak list
                    atomMass = CMass.get_atom_mass( obj.mMetal );
                    if max( data(:,1) ) < obj.mPrecursor - 0.02 + (atomMass - CMass.Proton)
                        aPeak = CPeak();
                        aPeak.mSpectrum = aSpectrum;
                        aPeak.mRawMZ = obj.mPrecursor + atomMass - CMass.Proton;
                        aPeak.mRawZ = 1;
                        aPeak.mIntensity = max(data(:,3)) + 10000;
                        if aSpectrum.mProtonated
                            aPeak.mMass = aPeak.mRawMZ;
                        end
                        aSpectrum.mPeaks(end+1) = aPeak;
                    end
                    
                    obj.mSpectrumSet(end+1) = aSpectrum;
                end
            end
        end
        
        function visualize_sources( obj )
            if isempty( obj.mEstimation ) || isempty( obj.mMassMap )
                return;
            end
            
            if obj.mSpectrumSet(1).mScanID > 0
                scanIdx = zeros(1, length(obj.mSpectrumSet));
                for k = 1 : length(scanIdx)
                    scanIdx(k) = obj.mSpectrumSet(k).mScanID;
                end
            else
                scanIdx = 1 : length(obj.mSpectrumSet);
            end
            
            colors = 'rbmgyk';
            x = 1 : size( obj.mMassMap, 1 );
            figure;

            num = length(obj.mSources.mean);
            legends = cell(1, num);
            for c = 1 : num
                legends{c} = ['Component ', num2str(c)];
            end
            
            cidx = 1;
            for c = 1 : num
                temp = exp(- (x-obj.mSources.mean(c)) .^ 2 / (2 * obj.mSources.std(c)^2) ) / obj.mSources.std(c);
                temp = temp .^ 2;
                plot(scanIdx, temp, colors(cidx) );
                if c == 1, hold on; end
                cidx = cidx + 1;
                if cidx > length(colors), cidx = length(colors); end
            end
            hold off;
            axis tight;
            title( 'Components' );
            legend( legends, 'FontSize', 11 );
        end % function visualize_sources
        
    end % method
    
    methods (Static)
        function [estPara, model] = mixture_decomposition_gaussian_kernel( massMap, numCom )
            % massmap[spectrumNumber, spectrumLength]
            % numCom is the number of components
            [spectrumNumber, spectrumLength] = size( massMap );
            xdata = (1 : spectrumNumber)';
            step = spectrumNumber/(numCom+1);
            start_point = zeros(1, 2*numCom + spectrumLength*numCom);
            lowBound = start_point;
            for k = 1 : numCom
                start_point(k*2-1) = k * step;
                start_point(k*2) = 2;
                lowBound(k*2-1) = 1;
                lowBound(k*2) = 2;
                sidx = floor((k-0.5)*step+1);
                eidx = floor((k+0.5)*step);
                for m = 1 : spectrumLength
                    start_point(2*numCom + (m-1)*numCom + k) = max( massMap(sidx:eidx, m) );
                end
            end
            model = @expfun;
            estPara = fmincon(model, start_point, [], [], [], [], lowBound );
            estPara( estPara < 0 ) = 0;

            function [sse, FittedCurve] = expfun(params)
                FittedCurve = zeros(spectrumNumber, spectrumLength);
                pp = zeros(spectrumNumber, numCom);
                for mm = 1 : numCom
                    % 1/\sigma * exp( -(x - \mu)^2 / (2 * \sigma^2) )
                    pp(:,mm) = exp(- (xdata-params(mm*2-1)) .^ 2 / (2 * params(mm*2)^2) ) / params(mm*2);
                end
                for kk = 1 : spectrumLength
                    offset = 2*numCom + (kk-1)*numCom;
                    FittedCurve(:,kk) = FittedCurve(:,kk) + pp * params(offset+1 : offset+numCom)';
                end
                ErrorVector = FittedCurve - massMap;
                paraPenalty = sum( params(2*numCom+1:end) > 0 );
                smoothPenalty = 0;
                offset = 2 * numCom;
                for kk = 2 : spectrumLength
                    dpara = params(offset+1 : offset+numCom) - params(offset+1+numCom : offset+numCom+numCom);
                    smoothPenalty = smoothPenalty + log(abs(dpara)+1);
                    offset = offset + numCom;
                end
                sse = sum(ErrorVector(:) .^ 2) + 5 * paraPenalty + sum(smoothPenalty);
            end
        end
        
        function [pcs, pcscores, recon, mu] = mixture_decomposition_pca( massmap, numCom )
            if nargin < 2 || isempty( numCom ), numCom = 2; end
            
            [pcs, pcscores, ~, ~, ~, mu] = pca( massmap' );
            pcs = pcs(:, 1:numCom);
            pcscores = pcscores(:, 1:numCom);
            recon = pcs * pcscores';
            recon = recon + mu' * ones(1, size(massmap, 2));
        end
        
        function vis_massmap( massmap, masses, transform )
            num = size( massmap, 2 );
            temp = sqrt(num / 12);
            rows = floor( temp * 3 );
            cols = ceil( num/rows );
            
            if nargin == 3 && ~isempty( transform )
                switch transform
                    case 'log'
                        massmap = log(massmap + 1);
                    case 'log2'
                        massmap = log2(massmap + 1);
                    case 'log10'
                        massmap = log10(massmap + 1);
                    case 'sqrt'
                        massmap = sqrt(massmap);
                end
            end
            
            maxx = size(massmap,1);
            maxy = max(max(massmap(:, 1:end-1)));
            figure;
            for k = 1 : num-1
                subplot(rows, cols, k);
                plot( massmap(:, k) );
                axis( [1, maxx, 0, maxy] );
                title( num2str(masses(k)) );
            end
            subplot(rows, cols, num);
            plot( massmap(:, num) );
            axis tight;
            title( num2str(masses(num)) );
        end
    end
end
