function str = vec2str(vec)
    if isempty(vec)
        str = '[]';
        return
    end
    str = '[';
    for k = 1 : (length(vec)-1)
        str = [str, num2str(vec(k)), ','];
    end
    str = [str, num2str(vec(end)), ']'];
end
