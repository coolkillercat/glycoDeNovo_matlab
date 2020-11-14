function [m, s, subSet] = robust_mean_std ( x, keep_percentage, verbose )
% function [m, s, subSet] = robust_mean_std ( x, keep_percentage, verbose )

if nargin < 2
    keep_percentage = 0.9;
end

if nargin < 3
    verbose = 0;
end

num = length( x );
% m = mean( x );
% s = std( x );
m = median( x );
s = sqrt( median( (x - m) .^ 2 ) );
preM = m - 10000;
preS = s - 10000;

preSubSet = x;
while ( abs( preM - m ) > abs( preM * 0.01 ) ) || ( abs( preS - s ) > preS * 0.01 )
    subSet = preSubSet( abs( preSubSet - m ) <= 3 * s );
    
    if length( subSet ) <= num * keep_percentage
        break;
    end

    preM = m;
    preS = s;
    %m = mean( subSet );
    %s = std( subSet );
    m = median( subSet );
    s = sqrt( median( (subSet - m) .^ 2 ) );

    if verbose
%       disp ( [preM, preS, m, s, length( preSubSet )] );
        disp ( [preM, preS, m, s]);
    end
    
    preSubSet = subSet;
end

m = mean( preSubSet );
s = std( preSubSet );