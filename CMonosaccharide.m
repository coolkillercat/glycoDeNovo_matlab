classdef CMonosaccharide < handle % By Pengyu Hong @ Brandeis University
    properties
        mID = -1;
        mSymbol = '';
        mName = '';
        mClass = '';
        mClassID = 0;
        mFormula = '';
        mInferredFormula = '';
        mConciseFormula = '';
        mPermethylated = 0;
        mReduced = 0;
        mMass = -1;
        
        mCarbon2Vertex = [];
        mDistance2Root = -1;
        mDistance2Leaf = -1;
        mVertices = struct( 'composition', {cell(1,6)}, 'mass', [0 0 0 0 0 0], 'modification', {cell(1,6)}, 'permethylated', [0 0 0 0 0 0] );
        
        mLegalLinkedInVs = [];
        mLinkedIn = struct( 'V', [], 'C', [], 'type', [], 'fromUnitID', [] ); % type: 1 - true; -1 - hypothetic.
        mLinkedTo = struct( 'V', [], 'C', [], 'type', [], 'toUnitID', [] ); % type: 1 - true; -1 - hypothetic.
        mCleavage = struct( 'type', '', 'side', '', 'vertices', [1 2 3 4 5 6], 'mass_to_mono', -1 );

        mInferredInLinkage = []; % struct( 'Vertex', [], 'Carbon', [], 'peaks', [], 'ions', {} ); % type: 1 - true; -1 - hypothetic.
        mInferredOutLinkage = []; % struct( 'Vertex', [], 'Carbon', [], 'peaks', [], 'ions', {} ); % type: 1 - true; -1 - hypothetic.
        mInferredInLinkageGroup = [] % struct( 'Vertex', [], 'Carbon', [], 'peaks', [] );
        mInferredOutLinkageGroup = [] % struct( 'Vertex', [], 'Carbon', [], 'peaks', [] );
    end
    
    methods
        function obj = CMonosaccharide( sym_name_class, permethylated )
            if nargin > 0
                if nargin < 2
                    permethylated = 0;
                end
                
                [ind, flag] = CMonosaccharideSet.find_monosaccharide_index( sym_name_class );
                if ~isempty( ind )
                    temp = CMonosaccharideSet.members(ind);
                    if flag
                        obj.mSymbol = temp.mSymbol;
                        obj.mName = temp.mName;
                    end
                    obj.mClass = temp.mClass;
                    obj.mClassID = temp.mClassID;
                    obj.mLegalLinkedInVs = temp.mLegalLinkedInVs;
                    obj.mCarbon2Vertex = temp.mCarbon2Vertex;
                    obj.mPermethylated = permethylated;
                    if obj.mPermethylated
                        obj.mMass = temp.mPermethylated.mass;
                        obj.mFormula = temp.mPermethylated.composition;
                        obj.mVertices.composition = temp.mPermethylated.vertex_composition;
                        obj.mVertices.mass = temp.mPermethylated.vertex_mass;
                        obj.mVertices.permethylated = temp.mPermethylated.vertex_permethylated;
                    else
                        obj.mMass = temp.mNative.mass;
                        obj.mFormula = temp.mNative.composition;
                        obj.mVertices.composition = temp.mNative.vertex_composition;
                        obj.mVertices.mass = temp.mNative.vertex_mass;
                    end
                else
                    throw( MException( 'CMonosaccharide:Constructor', ['Can not find the monossacharide: ', sym_name_class] ) );
                end
            end
        end
        
        function [linkedInCarbons] = add_linkedIn(obj, linkedInVs, linkedInTypes, fromUnitID )
            if isempty( linkedInVs )
                linkedInCarbons = [];
                return; 
            end
            if ~isempty( setdiff( linkedInVs, obj.mLegalLinkedInVs ) )
                throw( MException('CMonosaccharide:add_inV1', 'Try to link to an illegal vertex.') );
            end
            linkedInCarbons = zeros(1, length( linkedInVs ));
            for k = 1 : length( linkedInVs )
                if sum( obj.mLinkedIn.V == linkedInVs(k) ) > 0
                    throw( MException('CMonosaccharide:add_linkedIn', [num2str(linkedInVs(k)), '-th vertex occupied.']) );
                end
                temp = find( obj.mCarbon2Vertex == linkedInVs(k) );
                linkedInCarbons(k) = temp(end);
            end
            obj.mLinkedIn.V = [obj.mLinkedIn.V, linkedInVs];
            obj.mLinkedIn.C = [obj.mLinkedIn.C, linkedInCarbons];
            obj.mLinkedIn.type = [obj.mLinkedIn.type, linkedInTypes];
            obj.mLinkedIn.fromUnitID = [obj.mLinkedIn.fromUnitID, fromUnitID];
            if obj.mPermethylated
                obj.mVertices.modification( linkedInVs ) = {'-CH2'};
            end
            obj.update();
        end
        
        function add_inferred_inV(obj, linkedInVs, peaks)
            if nargin < 3, peaks = []; end
            changed = 0;
            if ~isempty( setdiff( linkedInVs, obj.mLegalLinkedInVs ) )
                throw( MException('CMonosaccharide:add_in_linkage', 'Try to link to an illegal vertex.') );
            end
            if isempty( obj.mInferredInLinkage )
                obj.mInferredInLinkage.V = linkedInVs;
                obj.mInferredInLinkage.peaks = peaks;
                changed = 1;
            else
                len = length( linkedInVs );
                for k = 1 : length( obj.mInferredInLinkage )
                    if sum( obj.mInferredInLinkage(k).V == linkedInVs ) == len
                        obj.mInferredInLinkage(k).peaks = unique( [obj.mInferredInLinkage(k).peaks, peaks] );
                        if length( obj.mInferredInLinkage(k).peaks ) > length( peaks )
                            changed = 1;
                        end
                        break;
                    end
                end
                if changed == 0
                    temp.V = linkedInVs; temp.peaks = peaks;
                    obj.mInferredInLinkage = [obj.mInferredInLinkage, temp];
                    changed = 1;
                end
            end
            if changed
                obj.update_inf_formula();
            end
        end
        
        function clear_linkedIn( obj )
            if ~isempty( obj.mLinkedIn.V )
                obj.mVertices.modification( obj.mLinkedIn.V ) = {[]};
                obj.mLinkedIn = struct( 'V', [], 'C', [], 'type', [], 'fromUnitID', [] );
                obj.update();
            end
        end
        function clear_linkedTo( obj )
            if ~isempty( obj.mLinkedTo.V )
                obj.mVertices.modification( obj.mLinkedTo.V ) = {[]};
                obj.mLinkedTo = struct( 'V', [], 'C', [], 'type', [], 'toUnitID', [] );
                obj.update();
            end
        end
        function remove_linkedIn( obj, V )
            idx = find( obj.mLinkedIn.V == V );
            if ~isempty( idx )
                flag = ones(1, length(obj.mLinkedIn.V));
                flag(idx) = 0; flag = flag > 0;
                obj.mLinkedIn.V = obj.mLinkedIn.V(flag);
                obj.mLinkedIn.C = obj.mLinkedIn.C(flag);
                obj.mLinkedIn.type = obj.mLinkedIn.type(flag);
            end
        end
        
        function [linkedInCarbons] = set_linkedIn (obj, linkedInVs, linkedInType, fromUnitID)
        % vtype = -1 if new inlinks are hypothetic
            if ~isempty( setdiff( linkedInVs, obj.mLegalLinkedInVs ) )
                throw( MException('CMonosaccharide:set_linkedIn', 'Try to link to the illegal carbon(s).') );
            end
            if nargin < 4
                fromUnitID = [];
            end
            
            % un-link
            unV = setdiff( obj.mLinkedIn.V, linkedInVs );
            if ~isempty( unV ) && obj.mPermethylated
                obj.mVertices.modification(unV) = {[]};
            end
            
            % link
            inV = setdiff( linkedInVs, obj.mLinkedIn.V );
            if ~isempty( inV ) && obj.mPermethylated
                obj.mVertices.modification(inV) = { '-CH2' };
            end
            
            obj.mLinkedIn.V = linkedInVs;
            obj.mLinkedIn.type = linkedInType;
            linkedInCarbons = zeros(1, length(linkedInVs));
            for b = 1 : length(linkedInVs)
                ids = find( obj.mCarbon2Vertex == linkedInVs(b) );
                linkedInCarbons(b) = ids(end);
            end
            obj.mLinkedIn.C = linkedInCarbons;
            obj.mLinkedIn.fromUnitID = fromUnitID;
            
            if obj.mPermethylated && ( ~isempty( unV )||~isempty( inV ) )
                obj.update();
            end
        end
        
        function set_linkedTo (obj, linkedToV, linkedToC, linkedToType, toUnitID )
            if obj.mID == 1
                throw( MException( 'CMonosaccharide:set_linkedTo', 'No linkout from the root!' ) );
            end
            obj.mLinkedTo.V = linkedToV;
            obj.mLinkedTo.C = linkedToC;
            obj.mLinkedTo.type = linkedToType;
            obj.mLinkedTo.toUnitID = toUnitID;
            if obj.mPermethylated && ~isempty(linkedToV)
                obj.mVertices.modification{2} = '-CH2'; % can only linked out from the 2nd vertex
                obj.update();
            end
        end
        
        function [x, a] = cleave( obj, crc_type )
            idx = find( strcmp( CCrossRingCleavage.cCRC, crc_type ) );
            if isempty( idx )
                throw( MException( 'CMonosaccharide:cleave', ['Wrong CRC: ', crc] ) );
            end
            a = obj.copy; 
            a.mCleavage.type = crc_type;
            a.mCleavage.side= 'A';
            a.mCleavage.vertices = CCrossRingCleavage.cVertex_A{idx};
            a.update();

            x = obj.copy; 
            x.mCleavage.type = crc_type;
            x.mCleavage.side = 'X';
            x.mCleavage.vertices = CCrossRingCleavage.cVertex_X{idx};
            x.update();
            
            a.mCleavage.mMass_to_mono = x.mMass;
            x.mCleavage.mMass_to_mono = a.mMass;
        end
        
        function complete( obj ) % If obj is cleaved, this function will add its complement.
            obj.mCleavage.type = '';
            obj.mCleavage.side = '';
            obj.mCleavage.mVertices = [1 2 3 4 5 6];
            obj.mCleavage.mMass_to_mono = -1;
            if obj.mPermethylated && obj.mID == 1
                obj.mVertices.modification{2} = 'CH2';
            end
            obj.update();
        end
        
        function mono = copy( obj )
            if isempty(obj)
                mono = [];
                return;
            end
            num = length(obj);
            mono = CMonosaccharide.empty( 0, num );
            for k = 1 : num
                mono(k) = CMonosaccharide;
                fields = fieldnames( obj(k) );
                for f = 1 : length( fields )
                    mono(k).(fields{f}) = obj(k).(fields{f});
                end
                mono(k).update();
            end
        end
        
        function set_reduced( obj, flag )
            if obj.mID == 1
                if flag ~= obj.mReduced
                    obj.mReduced = flag;
                    if obj.mReduced
                        if obj.mPermethylated
                            obj.mVertices.modification{1} = 'CH3';
                            obj.mVertices.modification{2} = 'H';
                        else
                            obj.mVertices.modification{1} = 'H';
                            obj.mVertices.modification{2} = 'H';
                        end
                    else
                        obj.mVertices.modification{1} = [];
                        obj.mVertices.modification{2} = [];
                    end
                    obj.update();
                end
            end
        end
        
        function set_permethylated( obj, flag )
            if obj.mPermethylated ~= flag
                obj.mPermethylated = flag;

                ind = [];
                if ~isempty( obj.mSymbol )
                    [ind, flag] = CMonosaccharideSet.find_monosaccharide_index( obj.mSymbol );
                elseif ~isempty( obj.mName )
                    [ind, flag] = CMonosaccharideSet.find_monosaccharide_index( obj.mName );
                elseif ~isempty( obj.mClass )
                    [ind, flag] = CMonosaccharideSet.find_monosaccharide_index( obj.mClass );
                end
                    
                if ~isempty( ind )
                    temp = CMonosaccharideSet.members(ind);
                    if flag
                        obj.mSymbol = temp.mSymbol;
                        obj.mName = temp.mName;
                    end
                    obj.mClass = temp.mClass;
                    if obj.mPermethylated
                        obj.mMass = temp.mPermethylated.mMass;
                        obj.composition = temp.mPermethylated.composition;
                        obj.mVertices.composition = temp.mPermethylated.vertex_composition;
                        obj.mVertices.mMass = temp.mPermethylated.vertex_mass;
                        obj.mVertices.modification = cell(1,6);
                        obj.mVertices.modification(obj.mLinkedTo.V) = '-CH2';
                        obj.mVertices.modification(obj.mLinkedIn.V) = {'-CH2'};
                        if obj.mID == 1  % is the root
                            obj.mVertices.modification{1} = 'CH3';
                            obj.mVertices.modification{2} = 'H';
                        end
                    else
                        obj.mMass = temp.native.mMass;
                        obj.composition = temp.native.composition;
                        obj.mVertices.composition = temp.native.vertex_composition;
                        obj.mVertices.mMass = temp.native.vertex_mass;
                        obj.mVertices.modification = cell(1,6);
                    end
                else
                    throw( MException( 'CMonosaccharide:Constructor', ['Can not find the monossacharide: ', obj.mSymbol] ) );
                end
            end
        end
            
        function set_symbol(obj, sym)
        	[ind, flag] = CMonosaccharideSet.find_monosaccharide_index( sym );
            if ~isempty( ind )
                temp = CMonosaccharideSet.members(ind);
                if flag
                    obj.mSymbol = temp.mSymbol;
                    obj.mName = temp.mName;
                end
                obj.mClass = temp.mClass;
                if obj.mPermethylated
                    obj.mMass = temp.mPermethylated.mMass;
                    obj.composition = temp.mPermethylated.composition;
                    obj.mVertices.composition = temp.mPermethylated.vertex_composition;
                    obj.mVertices.mMass = temp.mPermethylated.vertex_mass;
                else
                    obj.mMass = temp.native.mMass;
                    obj.composition = temp.native.composition;
                    obj.mVertices.composition = temp.native.vertex_composition;
                    obj.mVertices.mMass = temp.native.vertex_mass;
                end
            end
        end
        
        function update( obj )
            obj.mMass = sum( obj.mVertices.mass( obj.mCleavage.vertices ) );
            for k = 1 : length( obj.mCleavage.vertices )
                obj.mMass = obj.mMass + CMass.chem_formula_2_mass( obj.mVertices.modification{ obj.mCleavage.vertices(k) } );
            end
            obj.update_formula();
        end
           
        function update_formula( obj )
            obj.mFormula = '';
            obj.mInferredFormula = '';
            obj.mConciseFormula = '';
            
            if ~isempty( obj.mCleavage.type )
                obj.mFormula = [[obj.mCleavage.type, obj.mCleavage.side], '*'];
                obj.mInferredFormula = obj.mFormula;
                obj.mConciseFormula = obj.mFormula;
            end
            
            if ~isempty( obj.mSymbol )
                obj.mFormula = [obj.mFormula, obj.mSymbol];
                obj.mConciseFormula = [obj.mConciseFormula, obj.mSymbol];
            else
                obj.mFormula = [obj.mFormula, obj.mClass];
                obj.mConciseFormula = [obj.mConciseFormula, obj.mClass];
            end
            obj.mInferredFormula = [obj.mInferredFormula, obj.mClass];
        
            if ~isempty( obj.mLinkedTo ) && ~isempty( obj.mLinkedTo.type )
                if obj.mLinkedTo.C > 0
                    if obj.mLinkedTo.type == 1
                        obj.mFormula = [obj.mFormula, '-', num2str(obj.mLinkedTo.C)];
                    else
                        obj.mFormula = [obj.mFormula, '~', num2str(obj.mLinkedTo.C)];
                    end
                end
            end
            
            if ~isempty( obj.mInferredOutLinkage )
                temp = '';
                for k = 1 : length( obj.mInferredOutLinkage )
                    temp = [temp, sprintf( '%d(%d),', obj.mInferredOutLinkage(k).C, length( obj.mInferredOutLinkage(k).peaks ) )];
                end
                temp = temp(1:end-1);
                obj.mInferredFormula = [obj.mInferredFormula, '~', temp];
            end
        end
        
        function disp_inferred_linkage( obj )
            if ~isempty( obj.mInferredInLinkage )
                disp( 'Inferred inward linkage: ' );
                for k = 1 : length( obj.mInferredInLinkage )
                    cs = sprintf( '%d,', obj.mInferredInLinkage(k).C ); 
                    cs = cs(1:end-1);
                    disp( [obj.mClass, '<- ', cs, ' [', ...
                          num2str(length(obj.mInferredInLinkage(k).peaks)), ': ', ...
                          sprintf( '%d ', obj.mInferredInLinkage(k).peaks), ']' ] );
                end
            end
            if ~isempty( obj.mInferredOutLinkage )
                disp( 'Inferred outward linkage: ' );
                for k = 1 : length( obj.mInferredOutLinkage )
                    cs = sprintf( '%d,', obj.mInferredOutLinkage(k).C ); 
                    cs = cs(1:end-1);
                    disp( [obj.mClass, '-> ', cs, ' [', ...
                          num2str(length(obj.mInferredOutLinkage(k).peaks)), ': ', ...
                          sprintf( '%d ', obj.mInferredOutLinkage(k).peaks), ']' ] );
                end
            end
        end
        
        function ids = get_carbons( obj, vertexIDs )
            num = length(vertexIDs);
            ids = cell(1, num);
            for k = 1 : num
                ids{k} = find( obj.mCarbon2Vertex == vertexIDs(k) );
            end
        end
    end
    
    methods (Static)
        function link_monosaccharides( linkedToMono, linkedToV, linkedFromMono, linkType )
        % link the 2nd vertex of fromMono to linkedToV of toMono
        % linkType = -1 if hypothetic
            linkedToC = linkedToMono.add_linkedIn( linkedToV, linkType, linkedFromMono.mID );
            linkedFromMono.set_linkedTo( linkedToV, linkedToC, linkType, linkedToMono.mID );
        end
    end
end