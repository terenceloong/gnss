function plot_gnss_file(file_name)
% Plot first 0.1s data in specific file. (filename don't include suffix)

n = 4e5; %0.1s

[~,gnss_path] = system('echo %GNSS_PATH%');
if strcmp(gnss_path(1),'%')
    error('Can''t find environment variable GNSS_PATH !!!');
end
gnss_path(end) = '\';

file_path = [gnss_path, file_name, '.dat'];

fileID = fopen(file_path, 'r');
    data = fread(fileID, [2,n], 'int16'); %two row vector
    figure
    plot((1:n)/4e6, data(1,:))
    hold on
    plot((1:n)/4e6, data(2,:))
fclose(fileID);

end