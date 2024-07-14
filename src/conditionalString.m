function str = conditionalString(condition, trueStr, falseStr)
    if condition
        str = trueStr;
    else
        str = falseStr;
    end
end