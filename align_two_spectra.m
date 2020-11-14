function [matched, score] = align_two_spectra( spec1, spec2, massAccuracy, intensityStd )
% function [matched, cost] = align_two_spectra( spec1, spec2, massAccuracy, intensityStd )
% Align two spectra using Dynamic Programming

if nargin < 3
    massAccuracy = 0.01;
end
if nargin < 4
    intensityStd = 10000;
end

n1 = size( spec1, 1 );
n2 = size( spec2, 1 );
s = sparse( n1, n2 );
p = sparse( n1, n2 );

% Intialize the 1st row and 1st column.
for k = 1 : n2
    dMass = spec1(1,1) - spec2(k,1);
    if dMass > massAccuracy
        break;
    end
    if dMass < -massAccuracy
        continue;
    end   
    wm = 1 - abs(dMass);
    dIntensity = spec1(1,2) - spec2(k,2);
    wi = exp( -dIntensity / intensityStd );
    s(1,k) = wm * wi;
    p(1,k) = -1;
end
for k = 2 : n1
    dMass = spec1(k,1) - spec2(1,1);
    if dMass > massAccuracy
        break;
    end
    if dMass < -massAccuracy
        continue;
    end   
    wm = 1 - abs(dMass);
    dIntensity = spec1(k,2) - spec2(1,2);
    wi = exp( -dIntensity / intensityStd );
    s(k,1) = wm * wi;
    p(1,k) = -10;
end

for k = 2 : n1
    for m = 2 : n2
        dMass = spec1(k,1) - spec2(m,1);
        if dMass < -massAccuracy, continue; end
        if dMass > massAccuracy, break; end
        
        dIntensity = abs(spec1(k,2) - spec2(m,2));
        w = (1 - abs(dMass)) * exp(-dIntensity/intensityStd);
        
        o01 = s(k, m-1) + w - 0.5;
        o10 = s(k-1, m) + w - 0.5;
        o11 = s(k-1, m-1) + w;
        [~, idx] = max( [o01, o10, o11] );
        if idx == 1
            s(k,m) = o01; p(k,m) = -1;
        elseif idx == 2
            s(k,m) = o10; p(k,m) = -10;
        else
            s(k,m) = o11; p(k,m) = -11;
        end
    end
end

score = s(n1, n2);
matches