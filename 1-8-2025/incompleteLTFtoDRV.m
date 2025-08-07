% Open the source text file
inputFilename = 'Noise1to200Hz.ltf.txt';  % change as needed
fid = fopen(inputFilename, 'r');
if fid == -1
    error('Cannot open file: %s', inputFilename);
end

% Read and store header line
originalHeader = fgetl(fid);

% Read remaining lines as raw strings (to preserve formatting)
dataLines = {};
line = fgetl(fid);
while ischar(line)
    dataLines{end+1,1} = line*1e-3;
    line = fgetl(fid);
end
fclose(fid);

% Number of data entries
n = numel(dataLines);

% Generate time column starting at 0 with increment 0.005
time = (0:n-1)' * 0.005;

% Prepare output file and write header and data
outputFilename = 'Noise1to200Hz_convertable_0.txt';  % change as needed
fidOut = fopen(outputFilename, 'w');
if fidOut == -1
    error('Cannot open output file: %s', outputFilename);
end

% Write new header: time and the original header text
fprintf(fidOut, 'time \t PosT \t PosL \n');

% Write each line: time value and the original data line unchanged
for i = 1:n
    fprintf(fidOut, '%.6f\t %s\t 0.0 \n', time(i), dataLines{i});
end

fclose(fidOut);

fprintf('Wrote %d data lines to %s\n', n, outputFilename);
