s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';

hash = java.util.Hashtable;
C = {};
for i = 1:100
    key = randsample(s,20,true);
    value = randsample(s,40,true);
    hash.put(key, value);
    C{end+1} = {key, value};
end