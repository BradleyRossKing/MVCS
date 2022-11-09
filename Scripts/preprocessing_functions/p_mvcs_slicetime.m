function o_matlabbatch = p_mvcs_slicetime(files)

% Slice Time Correction
% BK; 2015/10/19

o_matlabbatch = [];

if size(files,1)>1 % Multiple files
    o_matlabbatch{end+1}.spm.temporal.st.scans = {cellstr(files)};
else % One 4D file
    o_matlabbatch{end+1}.spm.util.exp_frames.files = cellstr(files);
    o_matlabbatch{end}.spm.util.exp_frames.frames = Inf;
    dependancy = length(o_matlabbatch); % dependancy for next module
    o_matlabbatch{end+1}.spm.temporal.st.scans{1}(1) = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{dependancy}, '.','val', '{}',{dependancy}, '.','val', '{}',{dependancy}), substruct('.','files'));
end
o_matlabbatch{end}.spm.temporal.st.nslices = 43;
o_matlabbatch{end}.spm.temporal.st.tr = 2.65;
o_matlabbatch{end}.spm.temporal.st.ta = 2.5884; % TR-(TR/nslices)
o_matlabbatch{end}.spm.temporal.st.so = [1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42];
o_matlabbatch{end}.spm.temporal.st.refslice = 21;
o_matlabbatch{end}.spm.temporal.st.prefix = 'a_';

