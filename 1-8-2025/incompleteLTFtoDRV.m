% Open the source text file
inputFilename = 'Noise.ltf.txt';  % change as needed
fid = fopen(inputFilename, 'r');
if fid == -1
    error('Cannot open file: %s', inputFilename);
end

% Read and store header line
originalHeader = fgetl(fid);

% Read remaining lines as raw strings
dataLines = {};
line = fgetl(fid);
while ischar(line)
    % Store raw line (no arithmetic here)
    dataLines{end+1,1} = line;
    line = fgetl(fid);
end
fclose(fid);

% Number of data entries
n = numel(dataLines);

% Generate time column starting at 0 with increment 0.005
time = (0:n-1)' * 0.005;

% Prepare output file and write header and data
outputFilename = 'Noise_convertable_0.txt';  % change as needed
fidOut = fopen(outputFilename, 'w');
if fidOut == -1
    error('Cannot open output file: %s', outputFilename);
end

% Write new header: time and the original header text (plus any extra column names)
% Assuming originalHeader is single column name
fprintf(fidOut, 'time\t%s\tPosL\n', originalHeader);

% Write each line: time, scaled data value, and a zero constant (PosL)
for i = 1:n
    % Convert string to number and scale
    rawValue = str2double(dataLines{i});
    if isnan(rawValue)
        error('Non-numeric data encountered on line %d: "%s"', i+1, dataLines{i});
    end
    scaledValue = rawValue * 1e-3;
    fprintf(fidOut, '%.6f\t%.15g\t%.1f\n', time(i), scaledValue, 0.0);
end

fclose(fidOut);

fprintf('Wrote %d data lines to %s\n', n, outputFilename);
