function remove_git_folders_from_path()

c=textscan(path,'%s','delimiter',pathsep);
toberemoved = '';
for i=1:length(c{1,1})
    foo = cell2mat(c{1,1}(i));
    if regexp(foo,'\<\.git\>')
        toberemoved = cat(2,toberemoved,pathsep,foo);
    end
end
rmpath(toberemoved)
savepath