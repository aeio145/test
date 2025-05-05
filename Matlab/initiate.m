%Clear previous run:
clear; clc; close all;
parametersfile = 'parameters.txt';
valuesfile = 'values.csv';
datafile = 'data.csv';
% Check if parameters.csv file exists, if it doesn't stop the execution
if exist(parametersfile, 'file') == 2
    fprintf("'%s' file found\n", parametersfile);
else
    error("Error: The file '%s' was not found.", parametersfile);
end

% Check if values.csv file exists, if it doesn't stop the execution
if exist(valuesfile, 'file') == 2
    fprintf("'%s' file found\n", valuesfile);
else
    error("Error: The file '%s' was not found.", valuesfile);
end

% Check if data.csv file exists, if it doesn't stop the execution
if exist(datafile, 'file') == 2
    fprintf("'%s' file found\n", datafile);
else
    error("Error: The file '%s' was not found.", datafile);
end


% Check if Results folder exists. If it doesn't, create it.
if ~exist('Results', 'dir')
    mkdir('Results');
end

disp(">> Directory checks completed.");
fprintf('\n');
fprintf('\n');