c=textscan(path,'%s','delimiter',pathsep);
for i=1:length(c{1,1})
    foo = cell2mat(c{1,1}(i));
    if strfind(foo,'.svn')
        % fprintf('row %d is .svn folder!\n',i);
        rmpath(foo);
    end
end
savepath