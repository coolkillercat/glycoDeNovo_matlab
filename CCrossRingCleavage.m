classdef CCrossRingCleavage  % By Pengyu Hong @ Brandeis University
    
    % BA XY: B -> A + X <- Y
    properties (Constant)
        cCRC = {      '0,2',     '0,3',   '0,4',     '1,3',     '1,4',   '1,5',       '2,4',      '2,5',     '3,5' };
        cVertex_X = { [2 3],     [2 3 4], [2 3 4 5], [1 2 5 6], [1 2 6], [1 2],       [1 2 3 6],  [1 2 3],   [1 2 3 4] };
        cVertex_A = { [1 4 5 6], [1 5 6], [1 6],     [3 4],     [3 4 5], [3 4 5 6],   [4 5],      [4 5 6],   [5 6] };
    end
    
    methods (Static)
        function verify
            for k = 1 : length( CCrossRingCleavage.cCRC )
                if ~isempty( intersect( CCrossRingCleavage.cVertex_X{k}, CCrossRingCleavage.cVertex_A{k} ) )
                    disp( ['Vertex X n A: ', num2str(k)] );
                end
                if length( union( CCrossRingCleavage.cVertex_X{k}, CCrossRingCleavage.cVertex_A{k} ) ) ~= 6
                    disp( ['Vertex X u A: ', num2str(k)] );
                end
                if ~isempty( intersect( CCrossRingCleavage.cCarbon_X{k}, CCrossRingCleavage.cCarbon_A{k} ) )
                    disp( ['Carbon X n A: ', num2str(k)] );
                end
                if length( union( CCrossRingCleavage.cCarbon_X{k}, CCrossRingCleavage.cCarbon_A{k} ) ) ~= 6
                    disp( ['Carbon X u A: ', num2str(k)] );
                end
            end
        end
        
        function [x, a] = cleave( monosaccharide, crc_type )
        % function [x, a] = cleave( monosaccharide, crc_type )
            idx = find( strcmp( CCrossRingCleavage.cCRC, crc_type ) );
            if isempty( idx )
                throw( MException( 'CCrossRingCleavage:cleave', ['Wrong CRC: ', crc_type] ) );
            end
            a = monosaccharide.copy; 
            a.cleavage.type = crc_type; 
            a.cleavage.side = 'A';
            a.cleavage.vertices = CCrossRingCleavage.cVertex_A{idx};
            a.cleavage.mass = sum( a.vertices.mass( a.cleavage.vertices ) );

            x = monosaccharide.copy; 
            x.cleavage.type = crc_type; 
            x.cleavage.side = 'X';
            x.cleavage.vertices = CCrossRingCleavage.cVertex_X{idx};
            x.cleavage.mass = sum( x.vertices.mass( a.cleavage.vertices ) );
            
            a.cleavage.mass_to_mono = x.mass;
            x.cleavage.mass_to_mono = a.mass;
        end
        
        function info = get_CRC_info( crc_type )
            crcIdx = find(strcmp( CCrossRingCleavage.cCRCs, crc_type ) );
            if isempty(crcIdx)
                disp( ['Could not find the sub indexes for ', crc_type] );
            end
            info.crc_type = crc_type;
            info.vertex_X = CCrossRingCleavage.cVertex_X{ crcIdx };
            info.vertex_A = CCrossRingCleavage.cVertex_A{ crcIdx };
            info.carbon_X = CCrossRingCleavage.cCarbon_X{ crcIdx };
            info.carbon_A = CCrossRingCleavage.cCarbon_A{ crcIdx };
        end
        
        function monos = match_A( mass, err, preMe )
            if nargin < 2, err = 0.01; end
            if nargin < 3, preMe = 0; end
            monos = [];
            for k = 1 : length( CMonosaccharideSet.members )
                if preMe
                    for c = 1 : length( CCrossRingCleavage.cCRC )
                        temp = sum( CMonosaccharideSet.members(k).permethylated.vertex_mass( CCrossRingCleavage.cVertex_A ) );
                        if abs(temp - mass) < err
                            monos(end+1) = CMonosaccharide( CMonosaccharideSet.members(k).symbol, preMe );
                            monos(end).cleavage.type = CCrossRingCleavage.cCRC{c};
                            monos(end).cleavage.side = 'A';
                            monos(end).cleavage.vertices = CCrossRingCleavage.cVertex_A;
                            monos(end).cleavage.mass_to_mono = sum( CMonosaccharideSet.members(k).permethylated.vertex_mass( CCrossRingCleavage.cVertex_X ) );
                        end
                    end
                else
                end
            end
        end
        
        function monos = match_X( mass, err, preMe )
            if nargin < 2, err = 0.01; end
            if nargin < 3, preMe = 0; end
            monos = [];
            for k = 1 : length( CMonosaccharideSet.members )
                if preMe
                    for c = 1 : length( CCrossRingCleavage.cCRC )
                        temp = sum( CMonosaccharideSet.members(k).permethylated.vertex_mass( CCrossRingCleavage.cVertex_X ) );
                        if abs(temp - mass) < err
                            monos(end+1) = CMonosaccharide( CMonosaccharideSet.members(k).symbol, preMe );
                            monos(end).cleavage.type = CCrossRingCleavage.cCRC{c};
                            monos(end).cleavage.side = 'X';
                            monos(end).cleavage.vertices = CCrossRingCleavage.cVertex_X;
                            monos(end).cleavage.mass_to_mono = sum( CMonosaccharideSet.members(k).permethylated.vertex_mass( CCrossRingCleavage.cVertex_A ) );
                        end
                    end
                else
                end
            end
        end
    end
end