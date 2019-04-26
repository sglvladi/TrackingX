function output = isemptystruct(instruct)
%ISEMPTYSTRUCT Check if a given struct is empty (i.e. contains no fields)
    output = isempty(fieldnames(instruct));
end

