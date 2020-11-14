function new_mass = merge_masses( old_mass, threshold, advanced )

if nargin < 3
    advanced = 0;
end

tree = linkage( old_mass, 'complete' );

if advanced == 0
    c = cluster( tree, 'cutoff', threshold, 'criterion', 'distance' );
    newLen = max(c);
    new_mass = zeros(newLen, 1);
    
    for k = 1 : newLen
        tempM = old_mass( c == k );
        new_mass(k) = mean(tempM);
    end
    new_mass = sort(new_mass);
else
    n = length(old_mass);
    temp = cell(1, n*2 - 1);
    flag = ones(1, n*2 - 1);
    for k = 1 : n
        temp{k}.idxes = k;
        temp{k}.max_mass = old_mass(k);
        temp{k}.min_mass = old_mass(k);
        temp{k}.ave_mass = old_mass(k);
        temp{k}.diverse = 0;
    end
    for k = 1 : size(tree, 1)
        newIdx = n+k;
        temp{newIdx}.idxes = [temp{tree(k,1)}.idxes, temp{tree(k,2)}.idxes];
        temp{newIdx}.ave_mass = mean( old_mass(temp{newIdx}.idxes) );
        temp{newIdx}.max_mass = max( temp{tree(k,1)}.max_mass, temp{tree(k,2)}.max_mass );
        temp{newIdx}.min_mass = min( temp{tree(k,1)}.min_mass, temp{tree(k,2)}.min_mass );
        temp{newIdx}.diverse = max( [temp{newIdx}.max_mass, temp{newIdx}.min_mass] - temp{newIdx}.ave_mass );
        if temp{newIdx}.diverse >= threshold
            len = newIdx - 1;
            break;
        end
        flag(tree(k,1)) = 0;
        flag(tree(k,2)) = 0;
    end
    flag(len+1:end) = 0;
    temp = temp( flag > 0 );
    
    new_mass = zeros(1, length(temp));
    for k = 1 : length(temp)
        new_mass(k) = temp{k}.ave_mass;
    end
    new_mass = sort(new_mass);
end