function o_matlabbatch = p_coregister_multfiles(i_ref, i_source, i_2convert)
%
%   function o_matlabbatch = coregister(i_ref, i_source, i_2convert)
% 
%   i_ref:          [string]    Reference for coregistration  (steady)
%   i_source:       [string]    Source for coregistration
%   i_2coregister:  [string]    File to coregister
% 
%   abore: 17 Septembre 2015
%       - creation of coregister


o_matlabbatch = [];

if size(i_2convert,1)==1
    o_matlabbatch{end+1}.spm.util.exp_frames.files = cellstr(i_2convert);
    o_matlabbatch{end}.spm.util.exp_frames.frames = inf;
    o_matlabbatch{end}.spm.spatial.coreg.estimate.other = cfg_dep('Expand image frames: Expanded filename list.', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
else
    o_matlabbatch{end+1}.spm.spatial.coreg.estimate.other = cellstr(i_2convert); 
end
o_matlabbatch{end}.spm.spatial.coreg.estimate.ref = cellstr(i_ref);
o_matlabbatch{end}.spm.spatial.coreg.estimate.source = cellstr(i_source);
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 ...
0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
o_matlabbatch{end}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];
