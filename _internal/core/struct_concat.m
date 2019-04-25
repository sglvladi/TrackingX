function new_struct = struct_concat(struct1,struct2)
%STRUCT_CONCAT Concatenate two structures.

    f = fieldnames(struct2);
    for i = 1:length(f)
        struct1.(f{i}) = struct2.(f{i});
    end

    new_struct = struct1;

end

