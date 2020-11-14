classdef CMonosaccharideSet  % By Pengyu Hong @ Brandeis University
    properties (Constant)
        cMonoClasses = {'Xyl', 'Fuc', 'Hex', 'HexA', 'HexNAc', 'Kdo', 'Neu5Ac', 'Neu5Gc'};
        cNumberMonosaccharideClasses = 8;
    end
    
    methods (Static)
        function verify
            for k = 1 : length( CMonosaccharideSet.members )
                for m = 1 : length( CMonosaccharideSet.members(k).native.vertex_composition )
                    if abs( CMass.chem_formula_2_mass( CMonosaccharideSet.members(k).native.vertex_composition{m} ) - ...
                            CMonosaccharideSet.members(k).native.vertex_mass(m) ) > 0.001
                        disp( ['Native vertex mass does not match: ', CMonosaccharideSet.members(k).name, ' v-', num2str(m)] );
                    end
                end
                if abs( sum( CMonosaccharideSet.members(k).native.vertex_mass ) - ...
                             CMonosaccharideSet.members(k).native.mass) > 0.001
                    disp( ['Native mass does not match: ', CMonosaccharideSet.members(k).name] );
                end
                for m = 1 : length( CMonosaccharideSet.members(k).permethylated.vertex_composition )
                    if abs( CMass.chem_formula_2_mass( CMonosaccharideSet.members(k).permethylated.vertex_composition{m} ) - ...
                            CMonosaccharideSet.members(k).permethylated.vertex_mass(m) ) > 0.001
                        disp( ['Permethylated vertex mass does not match: ', CMonosaccharideSet.members(k).name, ' v-', num2str(m)] );
                    end
                end
                if abs( sum( CMonosaccharideSet.members(k).permethylated.vertex_mass ) - ...
                        CMonosaccharideSet.members(k).permethylated.mass) > 0.001
                    disp( ['Permethylated mass does not match: ', CMonosaccharideSet.members(k).name] );
                end
            end
        end
        
        function objs = get_mono_by_mass( mass, err, permethylated, byClass, possibleMonos )
            if nargin < 2 || isempty( err ), err = 0.01; end
            if nargin < 3 || isempty( permethylated ), permethylated = 0; end
            if nargin < 4 || isempty( byClass ), byClass = 1; end
            
            objs = []; 
            temp = [CMonosaccharideSet.members];
            if permethylated
                temp = [temp.permethylated];
            else
                temp = [temp.native];
            end
            temp = [temp.mass];
            ind = find( abs(temp-mass) <= err );
            if ~isempty( ind )
                objs = CMonosaccharide.empty(0, length(ind));
                for k = 1 : length(ind)
                    objs(k) = CMonosaccharide( CMonosaccharideSet.members(ind(k)).symbol, permethylated );
                end
            end
            
            if byClass
                num = length( objs );
                flag = ones(1, num);
                for a = 1 : num-1
                    for b = a+1 : num
                        if flag(b) && strcmp( objs(a).class, objs(b).class )
                            flag(b) = 0;
                        end
                    end
                end
                objs = objs(flag > 0);
                for a = 1 : length( objs )
                    objs(a).symbol = '';
                    objs(a).name = '';
                end
            end
            
            if ~isempty( possibleMonos )
                flag = zeros(1, length(objs));
                for a = 1 : length( objs )
                    flag(a) = sum( strcmp( {possibleMonos.class}, objs(a).class ) );
                end
                objs = objs( flag > 0 );
            end
        end
        
        function ms = get_mass( sym, permethylated, vertices )
        % sym - the symbol of CMonosaccharide.
        % permethylated - permethylated mass.
        % vertices - the masses of vertices. If vertices == [], return the whole mass
            if nargin < 2 || isempty(permethylated), permethylated = 0; end
            if nargin < 3, vertices = []; end
            
            ms = 0;
            ind = CMonosaccharideSet.find_monosaccharide_index( sym );
            
            if ~isempty(ind)
                if isempty( vertices )
                    if permethylated
                        ms = CMonosaccharideSet.members(ind).mPermethylated.mass;
                    else
                        ms = CMonosaccharideSet.members(ind).mNative.mass;
                    end
                else
                    if permethylated
                        ms = sum( CMonosaccharideSet.members(ind).mPermethylated.vertex_mass(vertices) );
                    else
                        ms = sum( CMonosaccharideSet.members(ind).mNative.vertex_mass(vertices) );
                    end
                end
            end
        end
        
        function inVs = get_allowed_inVs( sym )
            inVs = [];
            ind = CMonosaccharideSet.find_monosaccharide_index( sym );
            if ~isempty( ind )
                inVs = CMonosaccharideSet.members(ind).legalLinkedInVs;
            end
        end
        
        function [ind, flag] = find_monosaccharide_index( sym )
            temp = CMonosaccharideSet.members;
            flag = 1;
            ind = find( strcmp( {temp.mSymbol}, sym ) );
            if isempty( ind )
                ind = find( strcmp( {temp.mName}, sym ) );
                if isempty( ind )
                    flag = 0;
                    ind = find( strcmp( {temp.mClass}, sym ) );
                end
                if length( ind ) > 1
                    ind = ind(1);
                end
            end
        end
        
        function IDs = find_monosaccharide_classID( syms )
            IDs = zeros(1, length(syms));
            temp = CMonosaccharideSet.members;
            for k = 1 : length(syms)
                ind = find( strcmp( {temp.mClass}, syms{k} ) );
                if length( ind ) >= 1
                    ind = ind(1);
                    IDs(k) = CMonosaccharideSet.members(ind).mClassID;
                else
                    ind = find( strcmp( {temp.mSymbol}, syms{k} ) );
                    if length( ind ) >= 1
                        ind = ind(1);
                        IDs(k) = CMonosaccharideSet.members(ind).mClassID;
                    else
                        ind = find( strcmp( {temp.mName}, syms{k} ) );
                        if length( ind ) >= 1
                            ind = ind(1);
                            IDs(k) = CMonosaccharideSet.members(ind).mClassID;
                        end
                    end
                end
            end
        end
        
        function names = find_className_by_classID( IDs )
            temp = [CMonosaccharideSet.members.mClassID];
            names = cell(1, length(IDs));
            for k = 1 : length(IDs)
                idx = find(temp == IDs(k));
                if ~isempty(idx)
                    names{k} = CMonosaccharideSet.members(idx(1)).mClass;
                end
            end
        end
        
        function num = size()
            num = length( CMonosaccharideSet.members );
        end
    end
    
    properties (Constant)
        members = [struct( 'mSymbol', 'Xyl', ...
                           'mName', 'Xylose', ...
                           'mClass', 'Xyl', ...
                           'mClassID', 1, ...
                           'mLegalLinkedInVs', [3 4 5], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6], ...
                           'mNative', struct( 'mass', 150.0528234315, ...
                                             'composition', 'C5-H10-O5', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'CH2' }}, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 14.0156500642] ), ...
                           'mPermethylated', struct( 'mass', 206.1154236883, ...
                                                    'composition', 'C8-H16-O5', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'CH2' }}, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 0], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 14.0156500642] ) ), ...
                   struct( 'mSymbol', 'dHex', ...
                           'mName', 'Fucose', ...
                           'mClass', 'Fuc', ...
                           'mClassID', 2, ...
                           'mLegalLinkedInVs', [3 4 5], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 164.0684734957, ...
                                             'composition', 'C6-H12-O5', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H4'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 28.0313001284] ), ...
                           'mPermethylated', struct( 'mass', 220.1310737525, ...
                                                    'composition', 'C9-H18-O5', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 0], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 28.0313001284] ) ), ...
                    struct( 'mSymbol', 'Glc', ...
                           'mName', 'Glucose', ...
                           'mClass', 'Hex', ...
                           'mClassID', 3, ...
                           'mLegalLinkedInVs', [3, 4, 5, 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 180.0633881178, ...
                                             'composition', 'C6-H12-O6', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H4O'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 44.0262147505] ), ...
                           'mPermethylated', struct( 'mass', 250.1416384388, ...
                                                    'composition', 'C10-H20-O6', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C3H6O'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 58.0418648147] ) ), ...
                    struct( 'mSymbol', 'Gal', ...
                           'mName', 'Galactose', ...
                           'mClass', 'Hex', ...
                           'mClassID', 3, ...
                           'mLegalLinkedInVs', [3, 4, 5, 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 180.0633881178, ...
                                             'composition', 'C6-H12-O6', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H4O'} }, ... 
                                                                 'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 44.0262147505] ), ...
                           'mPermethylated', struct( 'mass', 250.1416384388, ...
                                                    'composition', 'C10-H20-O6', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C3H6O'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 58.0418648147] ) ), ...
                    struct( 'mSymbol', 'Man', ...
                           'mName', 'Mannose', ...
                           'mClass', 'Hex', ...
                           'mClassID', 3, ...
                           'mLegalLinkedInVs', [3, 4, 5, 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 180.0633881178, ...
                                             'composition', 'C6-H12-O6', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H4O'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 44.0262147505] ), ...
                           'mPermethylated', struct( 'mass', 250.1416384388, ...
                                                    'composition', 'C10-H20-O6', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C3H6O'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 58.0418648147] ) ), ...
                    struct( 'mSymbol', 'GlcA', ...
                           'mName', 'Glucuronic acid', ...
                           'mClass', 'HexA', ...
                           'mClassID', 4, ...
                           'mLegalLinkedInVs', [3, 4, 5], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 194.0426526757, ...
                                             'composition', 'C6-H10-O7', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H2O2'} }, ... 
                                                                 'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 58.0054793084] ), ...
                           'mPermethylated', struct( 'mass', 264.1209029967, ...
                                                    'composition', 'C10-H18-O7', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C3H4O2'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 72.0211293726] ) ), ...
                    struct( 'mSymbol', 'IdoA', ...
                           'mName', 'Iduronic acid', ...
                           'mClass', 'HexA', ...
                           'mClassID', 4, ...
                           'mLegalLinkedInVs', [3, 4, 5], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 194.0426526757, ...
                                             'composition', 'C6-H10-O7', ...
                                             'vertex_composition', {{'O', 'CH2O', 'CH2O', 'CH2O', 'CH2O', 'C2H2O2'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 30.0105646863 30.0105646863 30.0105646863 58.0054793084] ), ...
                           'mPermethylated', struct( 'mass', 264.1209029967, ...
                                                    'composition', 'C10-H18-O7', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C2H4O', 'C2H4O', 'C2H4O', 'C3H4O2'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 44.0262147505 44.0262147505 44.0262147505 72.0211293726] ) ), ...
                    struct( 'mSymbol', 'GalNAc', ...
                           'mName', 'N-acetylgalactosamine', ...
                           'mClass', 'HexNAc', ...
                           'mClassID', 5, ...
                           'mLegalLinkedInVs', [4 5 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 221.0899372193, ...
                                             'composition', 'C8-H15-N-O6', ...
                                             'vertex_composition', {{'O', 'CH2O', 'C3H5NO', 'CH2O', 'CH2O' ,'C2H4O'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 71.0371137878 30.0105646863 30.0105646863 44.0262147505] ), ...
                           'mPermethylated', struct( 'mass', 291.1681875403, ...
                                                    'composition', 'C12-H23-N-O6', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C4H7NO', 'C2H4O', 'C2H4O', 'C3H6O'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 85.0527638520 44.0262147505 44.0262147505 58.0418648147] ) ), ...
                    struct( 'mSymbol', 'GlcNAc', ...
                           'mName', 'N-acetylglucosamine', ...
                           'mClass', 'HexNAc', ...
                           'mClassID', 5, ...
                           'mLegalLinkedInVs', [4 5 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 3, 4, 5, 6, 6], ...
                           'mNative', struct( 'mass', 221.0899372193, ...
                                             'composition', 'C8-H15-N-O6', ...
                                             'vertex_composition', {{'O', 'CH2O', 'C3H5NO', 'CH2O', 'CH2O' ,'C2H4O'} }, ... 
                                             'vertex_mass', [15.9949146221 30.0105646863 71.0371137878 30.0105646863 30.0105646863 44.0262147505] ), ...
                           'mPermethylated', struct( 'mass', 291.1681875403, ...
                                                    'composition', 'C12-H23-N-O6', ...
                                                    'vertex_composition', {{'O', 'C2H4O', 'C4H7NO', 'C2H4O', 'C2H4O', 'C3H6O'} }, ... 
                                                    'vertex_permethylated', [0, 0, 1, 1, 1, 1], ...
                                                    'vertex_mass', [15.9949146221 44.0262147505 85.0527638520 44.0262147505 44.0262147505 58.0418648147] ) ), ...
                    struct( 'mSymbol', 'Kdo', ...
                           'mName', 'Kdo', ...
                           'mClass', 'Kdo', ...
                           'mClassID', 6, ...
                           'mLegalLinkedInVs', [4 5 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 2, 3, 4, 5, 6, 6, 6], ...
                           'mNative', struct( 'mass', 238.0688674262, ...
                                             'composition', 'C8-H14-O8', ...
                                             'vertex_composition', {{'O', 'C2H2O3', 'CH2', 'CH2O', 'CH2O', 'C3H6O2'} }, ... 
                                             'vertex_mass', [15.9949146221 74.0003939305 14.0156500642 30.0105646863 30.0105646863 74.0367794368] ), ...
                           'mPermethylated', struct( 'mass', 322.1627678114, ...
                                                    'composition', 'C14-H26-O8', ...
                                                    'vertex_composition', {{'O', 'C4H6O3', 'CH2', 'C2H4O', 'C2H4O', 'C5H10O2'} }, ... 
                                                    'vertex_permethylated', [0, 1, 0, 1, 1, 2], ...
                                                    'vertex_mass', [15.9949146221 102.0316940589 14.0156500642 44.0262147505 44.0262147505 102.0680795652] ) ), ...
                    struct( 'mSymbol', 'Neu5Ac', ...
                           'mName', 'N-acetylneuraminic acid', ...
                           'mClass', 'Neu5Ac', ...
                           'mClassID', 7, ...
                           'mLegalLinkedInVs', [4 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 2, 3, 4, 5, 6, 6, 6, 6], ...
                           'mNative', struct( 'mass', 309.1059812140, ...
                                             'composition', 'C11-H19-N-O9', ...
                                             'vertex_composition', {{'O', 'C2H2O3', 'CH2', 'CH2O', 'C3H5NO', 'C4H8O3'} }, ... 
                                             'vertex_mass', [15.9949146221 74.0003939305 14.0156500642 30.0105646863 71.0371137878 104.0473441231] ), ...
                           'mPermethylated', struct( 'mass', 407.2155316634, ...
                                                    'composition', 'C17-H31-N-O9', ...
                                                    'vertex_composition', {{'O', 'C4H6O3', 'CH2', 'C2H4O', 'C4H7NO', 'C7H14O3'} }, ... 
                                                    'vertex_permethylated', [0, 1, 0, 1, 1, 3], ...
                                                    'vertex_mass', [15.9949146221 102.0316940589 14.0156500642 44.0262147505 85.0527638520 146.0942943157] ) ), ...
                    struct( 'mSymbol', 'Neu5Gc', ...
                           'mName', 'N-glycolylneuraminic acid', ...
                           'mClass', 'Neu5Gc', ...
                           'mClassID', 8, ...
                           'mLegalLinkedInVs', [4 6], ... % allowed linked in vertices
                           'mCarbon2Vertex', [2, 2, 3, 4, 5, 6, 6, 6, 6], ...
                           'mNative', struct( 'mass', 325.1008958361, ...
                                             'composition', 'C11-H19-N-O10', ...
                                             'vertex_composition', {{'O', 'C2H2O3', 'CH2', 'CH2O', 'C3H5NO2', 'C4H8O3'} }, ... 
                                             'vertex_mass', [15.9949146221 74.0003939305 14.0156500642 30.0105646863 87.0320284099 104.0473441231] ), ...
                           'mPermethylated', struct( 'mass', 437.2260963497, ...
                                                    'composition', 'C17-H31-N-O10', ...
                                                    'vertex_composition', {{'O', 'C4H6O3', 'CH2', 'C2H4O', 'C5H9NO2', 'C7H14O3'} }, ... 
                                                    'vertex_permethylated', [0, 1, 0, 1, 1, 3], ...
                                                    'vertex_mass', [15.9949146221 102.0316940589 14.0156500642 44.0262147505 115.0633285383 146.0942943157] ) ), ...
                ];
    end
end