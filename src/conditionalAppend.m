function str = conditionalAppend(baseStr, condition, appendStr)
    if condition
        str = [baseStr appendStr];
    else
        str = baseStr;
    end
end