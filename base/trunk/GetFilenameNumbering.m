function [path fnprefix num_flag ext] = GetFilenameNumbering(filename)
% Return the file path, its name without ending numberings and its
% extension:
% [path fn num_flag ext]
% num_flag is the smallest number that the fnprefix is prefixed to.

[path fname ext]=fileparts(filename);
idx = regexpi([fname ext],['[0-9]+' ext '\>']);

fnprefix = fname;
if idx
    fnprefix = fname(1:(idx(end)-1));
    foo = dir([fullfile(path,fnprefix) '*' ext]);
    foo = foo(1).name;
    idx = regexpi(foo,['[0-9]+' ext '\>']);
    num_flag = str2num(foo(idx:end-4));
else
    num_flag=0;
end

