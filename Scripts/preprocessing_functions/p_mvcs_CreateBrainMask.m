function o_matlabbatch = p_mvcs_CreateBrainMask(segmented_anat, anat_folder)



o_matlabbatch = [];

o_matlabbatch{end+1}.spm.util.imcalc.input = cellstr(segmented_anat);
o_matlabbatch{end}.spm.util.imcalc.output = '';
o_matlabbatch{end}.spm.util.imcalc.outdir = cellstr(anat_folder);
o_matlabbatch{end}.spm.util.imcalc.expression = '(i1+i2+i3)';
o_matlabbatch{end}.spm.util.imcalc.var = struct('name', {}, 'value', {});
o_matlabbatch{end}.spm.util.imcalc.options.dmtx = 0;
o_matlabbatch{end}.spm.util.imcalc.options.mask = 0;
o_matlabbatch{end}.spm.util.imcalc.options.interp = 1;
o_matlabbatch{end}.spm.util.imcalc.options.dtype = 4;
o_matlabbatch{end+1}.spm.spatial.smooth.data(1) = cfg_dep('Image Calculator: ImCalc Computed Image: ', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
o_matlabbatch{end}.spm.spatial.smooth.fwhm = [2 2 2];
o_matlabbatch{end}.spm.spatial.smooth.dtype = 0;
o_matlabbatch{end}.spm.spatial.smooth.im = 0;
o_matlabbatch{end}.spm.spatial.smooth.prefix = 's';
o_matlabbatch{end+1}.spm.util.imcalc.input(1) = cfg_dep('Smooth: Smoothed Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
o_matlabbatch{end}.spm.util.imcalc.output = 'brain_mask';
o_matlabbatch{end}.spm.util.imcalc.outdir = cellstr(anat_folder);
o_matlabbatch{end}.spm.util.imcalc.expression = 'i1>0.01';
o_matlabbatch{end}.spm.util.imcalc.var = struct('name', {}, 'value', {});
o_matlabbatch{end}.spm.util.imcalc.options.dmtx = 0;
o_matlabbatch{end}.spm.util.imcalc.options.mask = 0;
o_matlabbatch{end}.spm.util.imcalc.options.interp = 1;
o_matlabbatch{end}.spm.util.imcalc.options.dtype = 4;