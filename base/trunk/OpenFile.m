function [fid st]=OpenFile(fn,att)
% [fid st]=OpenFile(fn,att)
% Opens a file whose file name is 'fn' and attributes of openning is 'att'
% It returns 1 if openning was NOT successful
st=0;
[fid,message]=fopen(fn,att);
if ~isempty(message)
    disp([message ': ' fn]);
    st=1;
    return
end