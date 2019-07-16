function fs_niftiribbon2itk(outfile, fillWithCSF)

if notDefined('fillWithCSF'), fillWithCSF = false; end

%% Convert freesurfer label values to itkGray label values
% We want to convert
%   Left white:   2 => 3
%   Left gray:    3 => 5
%   Right white: 41 => 4
%   Right gray:  42 => 6
%   unlabeled:    0 => 0 (if fillWithCSF == 0) or 1 (if fillWithCSF == 1)          

% read in the nifti
ni = readFileNifti(outfile);

% check that we have the expected values in the ribbon file
vals = sort(unique(ni.data(:)));
if ~isequal(vals, [0 2 3 41 42]')
    warning('The values in the ribbon file - %s - do no match the expected values [2 3 41 42]. Proceeding anyway...') %#ok<WNTAG>
end

% map the replacement values
invals  = [3 2 41 42];
outvals = [5 3 4  6];
labels  = {'L Gray', 'L White', 'R White', 'R Gray'};

fprintf('\n\n****************\nConverting voxels....\n\n');
for ii = 1:4;
    inds = ni.data == invals(ii);
    ni.data(inds) = outvals(ii);
    fprintf('Number of %s voxels \t= %d\n', labels{ii}, sum(inds(:)));
end

if fillWithCSF, ni.data(ni.data == 0) = 1; end

% write out the nifti
writeFileNifti(ni)

% done.
return
