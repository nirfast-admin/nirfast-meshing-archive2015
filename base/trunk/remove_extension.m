function [fn extension]=remove_extension(fname)
% [fn extension]=remove_extension(fname)
% removes the extension from string 'fn'
% it returns the removed extension too.

fn1=fname;
k=findstr(fn1,'.');
if ~isempty(k)
    k=k(length(k));
    extension = fn1(k:length(fn1));
    fn=fn1(1:k-1);
else
    fn=fname;
    extension=[];
end