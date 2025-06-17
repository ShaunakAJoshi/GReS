% Define the root folder (adjust if GReS is located elsewhere)
gresRoot = pwd;  % or use the full path if needed

% List of folders to check
foldersToCheck = { ...
    fullfile(gresRoot, "Code"), ...
    fullfile(gresRoot, "ThirdPartyLibs"), ...
    fullfile(gresRoot, "Utilities") ...
};

% Check if any of the folders is already on the MATLAB path
needRestore = false;
for k = 1:numel(foldersToCheck)
    if contains(path, foldersToCheck{k})
        needRestore = true;
        break;
    end
end

% Restore path only if needed
if needRestore
    restoredefaultpath;
end

% Add GReS paths
addpath(genpath("Code"));
addpath(genpath("ThirdPartyLibs"));
addpath(genpath("Utilities"));
% addpath(genpath("Tests_Develop"));

clear