% Import the parameters.txt file
params_data = readtable(parametersfile, 'Delimiter', ',');  % Comma-delimited file

% Extract variable names and values
var_names = params_data.Properties.VariableNames;  % Get variable names
values = table2array(params_data);  % Convert the table to an array of values

% Assign values to variables dynamically in the base workspace
for i = 1:length(var_names)
    assignin('base', var_names{i}, values(1, i));  % Assign value to variable in base workspace
    fprintf('%s = %.2f\n', var_names{i}, values(1, i));  % Display variable name and value
end

% Clear temporary variables to clean up the workspace
clear params_data var_names values i;
disp('Parameters acquired.');
fprintf('\n');




% Import the values.csv file
values_data = readtable('values.csv');

% Extract variable names and values
var_names = values_data.Properties.VariableNames;  % Get variable names
values = table2array(values_data);  % Convert the table to an array of values

% Assign values to variables dynamically in the base workspace
for i = 1:length(var_names)
    assignin('base', var_names{i}, values(1, i));  % Assign value to variable in base workspace
    fprintf('%s = %.2f\n', var_names{i}, values(1, i));  % Display variable name and value
end

% Clear temporary variables to clean up the workspace
clear values_data var_names values i;
disp('Values acquired.');
fprintf('\n');


% Import the data.csv file
data = readtable('data.csv');

% Extract variable names and values
var_names = data.Properties.VariableNames;  % Get variable names
values = table2array(data);  % Convert the table to an array of values (each column is a vector)

% Assign each column as a vector to variables in the base workspace
for i = 1:length(var_names)
    assignin('base', var_names{i}, values(:, i));  % Assign column as a vector to workspace
    
    % Display only the first value of each vector
    fprintf('%s = [%.2f, ...]\n', var_names{i}, values(1, i));  % Show only the first value
end

% Clear temporary variables to clean up the workspace
clear data var_names values i;
disp('Data acquired.');
fprintf('\n');


disp('>> All files successfully acquired.');
fprintf('\n');
fprintf('\n');

