function LTF_to_TXT( filename_with_ext, input_folder, varargin)
% LTF_to_TXT Convert a .drv, .tgt or .acq file via Python and load its .txt.
%
%   LTF_to_TXT(filename_with_ext)
%       Uses hard-coded Python module folder, module name, and function name.
%       filename_with_ext must include extension, e.g. 'MyRun.drv', 'MyRun.tgt' or 'MyRun.acq'.
%
%   LTF_to_TXT(input_folder, filename_with_ext,'OutputFolder', outF, 'Verbose', tf)
%       Optionally specify:
%         'OutputFolder' : folder where Python writes the .txt (default = InputFolder)
%         'Verbose'      : true/false for printing status messages (default: true)
%
%   The Python module folder, module name, and function name are fixed in the code:
%     PythonModuleFolder = '<hard-coded path>';
%     ModuleName = 'LTF_to_TXT';
%     FunctionName = 'ltf_to_txt';
%
%   The function:
%     1. Parses filename_with_ext, checks extension is .drv, .tgt, or .acq.
%     2. Builds full input path using InputFolder.
%     3. Ensures the Python module folder is on py.sys.path.
%     4. Imports the hard-coded module and calls the hard-coded function.
%     5. Expects the Python function to return full path to the generated .txt.
%     6. Calls loadTXT on that .txt and returns its result.

    %-------- Hard-coded Python settings --------
    PythonModuleFolder = 'C:\Users\afons\OneDrive - Universidade de Lisboa\Controlo de Plataforma Sismica\uniaxial_table_model\Adapting_Driver_Signal\personal_python_packages';
    ModuleName         = 'LTF_to_TXT';
    FunctionName       = 'ltf_to_txt';
    %--------------------------------------------

    %--- Parse inputs ---
    p = inputParser;
    p.KeepUnmatched = false;
    p.addRequired('filename_with_ext', @(x) ischar(x) || isstring(x));
    p.addParameter('OutputFolder', '', @(x) ischar(x) || isstring(x));
    p.addParameter('Verbose',      true, @(x) islogical(x) || isnumeric(x));
    p.parse(filename_with_ext, varargin{:});

    filename_with_ext = char(p.Results.filename_with_ext);
    output_folder    = char(p.Results.OutputFolder);
    verbose          = logical(p.Results.Verbose);


    if isempty(output_folder)
        output_folder = input_folder;
        if verbose
            fprintf('Using default OutputFolder (same as InputFolder): %s\n', output_folder);
        end
    end

    %--- Parse filename and extension ---
    [~, base_name, ext] = fileparts(filename_with_ext);
    ext = lower(ext);
    if verbose
        fprintf('Parsed filename: base="%s", extension="%s"\n', base_name, ext);
    end

    %--- Accept .drv, .tgt, or .acq ---
    if ~ismember(ext, {'.drv', '.tgt', '.acq'})
        error('Unsupported extension "%s". Must be ".drv", ".tgt", or ".acq".', ext);
    end

    %--- Construct full input path and check existence ---
    in_file = fullfile(input_folder, filename_with_ext);
    if verbose
        fprintf('Constructed input file path: %s\n', in_file);
    end
    if ~isfile(in_file)
        error('Input file does not exist: %s', in_file);
    end

    %--- Ensure Python module folder is on py.sys.path ---
    try
        already = any(cellfun(@(p) strcmp(char(p), PythonModuleFolder), cell(py.sys.path)));
        if ~already
            insert(py.sys.path, int32(0), PythonModuleFolder);
            if verbose
                fprintf('Inserted Python module folder into sys.path: %s\n', PythonModuleFolder);
            end
        elseif verbose
            fprintf('Python module folder already on sys.path.\n');
        end
    catch ME
        id = ME.identifier;
        if isempty(id)
            id = 'MyApp:PythonPathError';
        end
        warning(id,'Could not modify Python sys.path: %s', ME.message);
    end

    %--- Import Python module ---
    try
        pyModule = py.importlib.import_module(ModuleName);
        if verbose
            fprintf('Imported Python module: %s\n', ModuleName);
        end
    catch ME
        fprintf('Error importing Python module "%s": %s\n', ModuleName, ME.message);
        if isprop(ME, 'PYTHONERROR')
            fprintf('Python traceback:\n%s\n', char(ME.PYTHONERROR));
        end
        error('Cannot proceed without the Python module.');
    end

    %--- Call the Python conversion function ---
    try
        py_output = pyModule.(FunctionName)(in_file, output_folder);
        if verbose
            fprintf('Called Python %s.%s on %s\n', ModuleName, FunctionName, in_file);
        end

        if isempty(py_output) || isequal(py_output, py.None)
            warning('Python function returned None or empty. Check if conversion succeeded in folder: %s', output_folder);
            output_path = '';
        else
            output_path = char(py_output);
            if verbose
                fprintf('Python returned output path: %s\n', output_path);
            end
        end
    catch ME
        fprintf('Error calling Python function %s.%s: %s\n', ModuleName, FunctionName, ME.message);
        if isprop(ME, 'PYTHONERROR')
            fprintf('Python traceback:\n%s\n', char(ME.PYTHONERROR));
        end
        error('Python conversion failed.');
    end

    %--- Ensure output_path is available ---
    if isempty(output_path)
        error('No output path from Python; cannot load TXT.');
    end
    if ~isfile(output_path)
        warning('Python reported output path but file not found: %s', output_path);
    end

end
