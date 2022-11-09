function o_matlabbatch = p_mvcs_3Dto4D(epi)



o_matlabbatch = [];

o_matlabbatch{end+1}.spm.util.cat.vols = cellstr(epi); 
o_matlabbatch{end}.spm.util.cat.name = '4D_ra_vol.nii';
o_matlabbatch{end}.spm.util.cat.dtype = 4;
o_matlabbatch{end}.spm.util.cat.RT = NaN;