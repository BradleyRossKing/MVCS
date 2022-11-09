function o_matlabbatch = p_mvcs_reverse_norm(flow, ROI)

o_matlabbatch = [];

o_matlabbatch{end+1}.spm.spatial.normalise.write.subj.def = cellstr(flow);
o_matlabbatch{end}.spm.spatial.normalise.write.subj.resample = cellstr(ROI);
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
                                                          78 76 85];
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.vox = [2 2 2];
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.interp = 4;
o_matlabbatch{end}.spm.spatial.normalise.write.woptions.prefix = 'w';