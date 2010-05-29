function systemcommand = GetSystemCommand(command_name)
% systemcommand = GetSystemCommand(command_name)
% Determines the platform that matlab is running and calls the appropriate
% executable file using following convention:
% command_name-mac.exe  For Mac
% command_name-linux.exe  For Linux
% command_name.exe For Windows
% 
% If command_name is not found an empty string will be returned.
% 
% Note that systemcommand my contain spaces so the calling functin should
% enclose it using double quotes.

os=computer;
systemcommand=[];
if ~isempty(strfind(os,'PCWIN')) % Windows
    systemcommand = which([command_name '.exe']);
elseif ~isempty(strfind(os,'MAC')) % Mac OS
    systemcommand = which([command_name '-mac.exe']);
elseif ~isempty(strfind(os,'GLNX')) % Linux
    systemcommand = which([command_name '-linux.exe']);
else
    fprintf('\n\tCan not establish your platform!\n');
    error('\t%s will not be executed\n\n', command_name);
end
if isempty(systemcommand)
    error('Could not find ''%s'' !', command_name);
end


