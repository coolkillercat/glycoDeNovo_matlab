%% Get all monosaccharides
monoClasses = unique( {CMonosaccharideSet.members.mClass} );
num = length( monoClasses );

%%
permethylation = 1;
monos = CMonosaccharide.empty(0, num);
for k = 1 : num
    monos(k) = CMonosaccharide( monoClasses{k}, permethylation );
end

%%
numCRC = length( CCrossRingCleavage.cCRC );
fieldCRC = cell(1, numCRC);
for k = 1 : numCRC
    % fieldCRC{k} = ['crc_', strrep(  CCrossRingCleavage.cCRC{k}, ',', '' )];
    fieldCRC{k} = strrep( CCrossRingCleavage.cCRC{k}, ',', '' );
end

%% collect A cross-rings
crcSet = [];
for k = 1 : num
    % crcSet.(monoClasses{k}) = [];
    for m = 1 : numCRC        
        [X, A] = monos(k).cleave( CCrossRingCleavage.cCRC{m} );
        crcSet.([monoClasses{k}, fieldCRC{m}, 'A']) = A;
        crcSet.([monoClasses{k}, fieldCRC{m}, 'X']) = X;
    end
end

%% mono + crc
mono_crc = [];
crcSetItems = fields( crcSet );
ind = 1;
for k = 1 : num
    for m = 1 : length(crcSetItems)
        if crcSetItems{m}(end) == 'X', continue; end
        mono_crc(ind).crc = [monos(k).mClass, ' + ', crcSetItems{m}];
        mono_crc(ind).topo = [monoClasses{k}, '_', crcSetItems{m}];
        mono_crc(ind).mass = monos(k).mMass + crcSet.(crcSetItems{m}).mMass - CMass.H2O - 3 * permethylation * CMass.CH2;
        ind = ind + 1;
    end
end

%% find common X CRCs
fprintf( 'CRC' );
fprintf( '\t%s', monoClasses );
for m = 1 : numCRC
    masses = zeros(1, num);
    for k = 1 : num
        masses(k) = crcSet.([monoClasses{k}, '_', fieldCRC{m}, '_A']).mMass;
    end
    fprintf( '%s', fieldCRC{m} );
    fprintf( '\t%6.3f', masses );
    fprintf( '\n' );
end

%%
if permethylation
    disp( '--- Permethylated B-ion: Monosaccharide == Cross-Ring Fragment ---' );
else
    disp( '--- B-ion: Monosaccharide == Cross-Ring Fragment ---' );
end

for k = 1 : num
    monoMass = monos(k).mMass - CMass.H2O - 2 * permethylation * CMass.CH2;
    for m = 1 : length( crcSetItems )
        if crcSetItems{m}(end) == 'X', continue; end
        diff = monoMass - (crcSet.(crcSetItems{m}).mMass - CMass.H2O - permethylation * CMass.CH2);
        if abs( diff ) < 0.005
            disp( [monos(k).mClass, ' == ', crcSet.(crcSetItems{m}).mFormula, ': mass diff = ', num2str(diff)] );
        end
    end
end
disp( '-------------------' );

%%
if permethylation
    disp( '--- Permethylated B-ion: Monosaccharide == Monosaccharide + Cross-Ring Fragment ---' );
else
    disp( '--- B-ion: Monosaccharide == Monosaccharide + Cross-Ring Fragment ---' );
end

mono_cr_masses = [mono_crc.mass];
for k = 1 : num
    diff = mono_cr_masses - (monos(k).mMass - CMass.H2O - 2 * permethylation * CMass.CH2);
    ind = find( abs(diff) < 0.005 );
    for m = 1 : length(ind)
        disp( [monos(k).mClass, ' == ', mono_crc(ind(m)).topo, ': ', num2str(diff(ind(m)))] );
    end
end
disp( '-------------------' );
