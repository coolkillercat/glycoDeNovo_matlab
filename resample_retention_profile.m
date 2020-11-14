function result = resample_retention_profile( data, osize )

result = zeros(1, osize);
p = data / sum(data);
sp = cumsum(p);
num = round( osize * p );
diff = osize - sum(num);

if diff ~= 0
    s = sign(diff);
    diff = abs(diff);
    k = 1;
    while k <= diff
        a = rand;
        [~, idx] = max( find(sp <= a) );
        if s < 0
            if num(idx) == 0
                continue;
            else
                num(idx) = num(idx) - 1;
            end
        else
            num(idx) = num(idx) + 1;
        end
        k = k + 1;
    end
end

result = [];
for k = 1 : length(num)
    result = [result, ones(1, num(k))*k];
end