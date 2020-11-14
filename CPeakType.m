classdef CPeakType
    properties (Constant)
        Unknown = 0;
        B = 1;
        C = bi2de([0 1]);
        A = bi2de([0 0 1]);
        Y = bi2de([0 0 0 1]);
        Z = bi2de([0 0 0 0 1]);
        X = bi2de([0 0 0 0 0 1]);
        T = bi2de([0 0 0 0 0 0 1]);
        Minus2H = bi2de([0 0 0 0 0 0 0 1]);
        
        B_H = 11;
        B_2H = 12;
        C_H = 21;
        C_2H = 22;
        T_H = 31;
        T_2H = 32;
    end
    
    methods(Static)
        function r = isA( intype )
            r = bitand(intype, CPeakType.A);
        end
        
        function r = isB( intype )
            r = bitand(intype, CPeakType.B);
        end
        
        function r = isC( intype )
            r = bitand(intype, CPeakType.C);
        end
        
        function r = isX( intype )
            r = bitand(intype, CPeakType.X);
        end
        
        function r = isY( intype )
            r = bitand(intype, CPeakType.Y);
        end
        
        function r = isZ( intype )
            r = bitand(intype, CPeakType.Z);
        end
        
        function r = isMinus2H( intype )
            r = bitand(intype, CPeakType.Minus2H);
        end
        
        function r = combineType( intype, addtype )
            r = bitor( intype, addtype );
        end
        
        function result = checkType( intype )
            result = '';
            if CPeakType.isA(intype)
                result = [result, 'A'];
            elseif CPeakType.isB(intype)
                result = [result, 'B'];
            elseif CPeakType.isC(intype)
                result = [result, 'C'];
            elseif CPeakType.isX(intype)
                result = [result, 'X'];
            elseif CPeakType.isY(intype)
                result = [result, 'Y'];
            elseif CPeakType.isZ(intype)
                result = [result, 'Z'];
            end
            
            if CPeakType.isMinus2H( intype )
                result = [result, '-2H'];
            end
        end
    end
end