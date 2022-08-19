
% NECESSITATES SPM TO RUN
% CALLS ROI FILES (IN NATIVE SPACE) BUILT IN SPM AND FSL
% CALLS IMAGE FILES PREPROCESSED IN SPM
% CALLS OUTPUT OF OTHER SCRIPT SCRUBBING PRECHECK


% DIRECTORIES OF INTEREST
org.main_folder = {'/Volumes/DATA/MVCS_ALLO_EGO_20200406/MVCS/DATA_5runs'};
org.analysesDate = {'_analyses_20211022'}; % output data / figures in folder with this Date


% LIST OF ROIs (must match name of .nii file stored in each individual' T1 folder)
org.roi = { ...
    'wAE_M1_Right_Intersect_FWE005', ...
    'wBil_Seg_Puta', ...
    'wBil_Seg_Hipp', ...
    }';


% PARAMETERS
org.identifier_func = 'r_a'; % prefix of .img containing functional data
org.identifier_gm = 'c1s'; % prefix of .nii with GM segmentation
org.identifier_wm = 'c2s'; % prefix of .nii with WM segmentation
org.identifier_csf = 'c3s'; % prefix of .nii with CSF segmentation
org.reference_run = 1; % functional run to serve as a reference for the affine transformation matrix; needs to be a run that is valid (present) for each individual; 
org.head_size=50;  % diameter of a standard head (in mm)
org.TR = [2.65; 2.65; 2.65; 2.65; 2.65]; % TR for each run
org.TR_ref = 2.65; % Setting the reference TR (resampled such that all data acquired with this TR)
org.preprocessing.num_of_pcs = 3; % # of PCs extracted from PCA on WM and CSF signal to use as nuisance regressors
org.preprocessing.erodeROI = [0 0 0]; % 0 if step not done; 1 if done; input for each ROI
org.preprocessing.scrubbing_fd_thres=0.5; % threshold (in mm) for determining bad volumes based on framewise displacement 
org.preprocessing.scrubbing = 1; % 0 if step not done; 1 if done
org.preprocessing.scrubbing_backwards = 0; % number of volumes to be removed BEFORE volume detected to be bad
org.preprocessing.scrubbing_forwards = 1;  % number of volumes to be removed AFTER volume detected to be bad
org.preprocessing.wm_regression = 1; % 0 if step not done; 1 if done
org.preprocessing.csf_regression = 1; % 0 if step not done; 1 if done
org.preprocessing.move_regression = 1; % 0 if step not done; 1 if done
org.preprocessing.gm_check = 1; % 0 if step not done; 1 if done; voxel within ROI must be have a minimum probability of being considered GM to be included (probability specified next)
org.preprocessing.gm_check_prob = 0.1; %  voxel within ROI must be have at least 0.10 prob of being considered GM to be included
org.preprocessing.volume_number_match = 1; % 0 if step not done; 1 if done; matches the number of volumes across all runs within a participant (after scrubbing and re-sampling for different TR)
org.preprocessing.hp_filter = 1; % 0 if step not done; 1 if done
org.preprocessing.hp_value = (1/128); 
load('colors_mvcs.mat');

% IF MATCHING # OF VOLUMES ACROSS RUN, LOAD THE .MAT FILE THAT SPECIFIES
% RESULTS OF A SCRUBBING PRECHECK
if org.preprocessing.volume_number_match == 1
    load([char(org.main_folder) '/_scrubbing_precheck/' 'scrubbing_precheck' char(org.analysesDate) '.mat'], 'precheck');
end

files = dir(char(org.main_folder));
counter = 0;
for ii = 1:1:length(files)
    if contains(files(ii).name,'ALLO_EGO')
        counter = counter + 1;
        org.subjects(counter,1) = {files(ii).name};    
    end
end
clear counter; clear ii; clear files;

org.sessions = repmat({'NaN'},length(org.subjects),5);  % 5 is chosen as a max 
for sub = 1:1:length(org.subjects)
    files = dir([char(org.main_folder) '/' char(org.subjects(sub))]);
    counter = 0;
    for ii = 1:1:length(files)
        if strcmp(files(ii).name,'.') || strcmp(files(ii).name,'..')
        elseif contains(files(ii).name,'ep2d')
            counter = counter + 1;
            org.sessions(sub,counter) = {files(ii).name};
        elseif contains(files(ii).name,'MPRAGE')
            org.anat(sub,1) = {files(ii).name};
        end
    end
    clear files; clear counter; clear ii;
end
clear sub;

%%

for nSub =1:length(org.subjects) % loop through subjects
    % ACCESS INFORMATION FROM A REPRESENTATIVE (REFERENCE) FUNCTIONAL RUN OF THIS SUBJECT
    ref_func_folder = fullfile(org.main_folder{1}, org.subjects{nSub},org.sessions{nSub,org.reference_run}, filesep); 
    ref_func = p_mvcs_get_files(ref_func_folder,'img', org.identifier_func); 
    ref_func = cellstr(ref_func);
    ref_func=ref_func{1};
    Vref=spm_vol(ref_func);
    ref_mat=Vref.mat; % get transform matrix mapping voxel coordinates to real world coordinates of a exemplar functional run
    ref_dim(1)=Vref.dim(1);
    ref_dim(2)=Vref.dim(2);
    ref_dim(3)=Vref.dim(3);
    A=ref_mat(1:3,1:3);
    B=ref_mat(1:3,4);
    
    clear brain_mask; clear brain_mask_index; 
    anat_folder = fullfile(org.main_folder{1}, org.subjects{nSub},org.anat{nSub},filesep);
    if not(isfolder([anat_folder char(org.analysesDate) '/']))
        mkdir([anat_folder char(org.analysesDate) '/'])
    end
    
    % Obtain Brain Mask
    disp(['brain mask for subject ' num2str(nSub)]);
    copyfile([anat_folder 'brain_mask.nii' ],[anat_folder char(org.analysesDate) '/' 'brain_mask.nii' ]);
    brain_mask = p_mvcs_get_files([anat_folder char(org.analysesDate) '/'],'nii','brain_mask'); % Find brain mask of subject (native space)
    Vt=spm_vol(brain_mask);
    mat_mask=Vt.mat; % get transform matrix mapping voxel coordinates to real world coordinates of the brain mask
    if sum(abs(ref_mat(:)-mat_mask(:)))>0 % if mismatch between affine transform matrices from mask and functional references; will reslice brain mask and over-write image    
        tmp=spm_read_vols(Vt);
        Vt.pinfo=[0.001;0;0];
        Vt.dt=[8 0];
        spm_write_vol(Vt,tmp>0);
        clear tmp; 

        clear matlabbatch;
        matlabbatch{1}.spm.spatial.coreg.write.ref{1} = ref_func;
        matlabbatch{1}.spm.spatial.coreg.write.source{1} = brain_mask; 
        matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
        matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
        matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
        spm_jobman('run', matlabbatch);
        clear matlabbatch; 
        
        Vtx=spm_vol([anat_folder char(org.analysesDate) '/' 'rbrain_mask.nii' ]);
        tmp=spm_read_vols(Vtx);
        Vtx.fname = brain_mask;
        delete([anat_folder char(org.analysesDate) '/brain_mask.nii']);
        delete([anat_folder char(org.analysesDate) '/rbrain_mask.nii']);
        spm_write_vol(Vtx,tmp);
        clear Vtx; clear tmp;
    end % if sum(abs(ref_mat(:)-mat_mask(:)))>0 
    clear Vt; clear mat_mask; 
    Vt=spm_vol(brain_mask);
    mask=spm_read_vols(Vt);
    brain_mask_index=find(reshape(mask,1,numel(mask)) > 0.5*max(mask(:)));
    nvoxel=length(brain_mask_index);
    clear Vt; clear mask; 
    
    % Identify x, y, z positions of voxels
    zpos=ceil(brain_mask_index/(ref_dim(1)*ref_dim(2)));
    ypos=ceil((brain_mask_index-(zpos-1)*ref_dim(1)*ref_dim(2))/ref_dim(1));
    xpos=brain_mask_index-(zpos-1)*ref_dim(1)*ref_dim(2)-(ypos-1)*ref_dim(1);
    voxel_position=A*[xpos; ypos; zpos]+B*ones(1,nvoxel);
    clear zpos; clear ypos; clear xpos; 
    
    % White Matter Mask - used to extract nuisance regressors
    if org.preprocessing.wm_regression == 1
        disp(['WM mask for subject ' num2str(nSub)]);
        wm_mask = p_mvcs_get_files(anat_folder,'nii', org.identifier_wm); % Find brain mask of subject (native space)
        copyfile(wm_mask,[anat_folder char(org.analysesDate) '/' 'wm_mask.nii' ]);
        wm_mask = p_mvcs_get_files([anat_folder char(org.analysesDate) '/'],'nii','wm_mask'); 
        Vt=spm_vol(wm_mask);
        tmp=spm_read_vols(Vt);
        Vt.pinfo=[0.001;0;0];
        Vt.dt=[8 0];

        if sum(abs(ref_mat(:)-Vt.mat(:)))>0 % if mismatch between affine transform matrices from mask and functional references; will reslice brain mask and over-write image    
            spm_write_vol(Vt,tmp);

            clear matlabbatch;
            matlabbatch{1}.spm.spatial.coreg.write.ref{1} = ref_func;
            matlabbatch{1}.spm.spatial.coreg.write.source{1} = Vt.fname; 
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            spm_jobman('run', matlabbatch);
            clear matlabbatch; 

            Vtx=spm_vol([anat_folder char(org.analysesDate) '/' 'rwm_mask.nii' ]);
            tmp=spm_read_vols(Vtx);
            Vtx.fname = wm_mask;
            delete([anat_folder char(org.analysesDate) '/wm_mask.nii']);
            delete([anat_folder char(org.analysesDate) '/rwm_mask.nii']);
            spm_write_vol(Vtx,tmp);
            clear Vtx;   

        end % if sum(abs(ref_mat(:)-Vt.mat(:)))>0   
        tmp=tmp>=prctile(tmp(tmp>0),90); % take those that are greater than 90th percentile (essentially those with WM prob = 1)
        mask=tmp(brain_mask_index);
        wm_mask_index=find(mask > 0.5*max(mask(:)));
        clear tmp; clear mask; clear Vt; 
    end
    
    % CSF Mask - used to extract nuisance regressors
    if org.preprocessing.csf_regression == 1
        disp(['CSF mask for subject ' num2str(nSub)]);
        csf_mask = p_mvcs_get_files(anat_folder,'nii', org.identifier_csf); % Find brain mask of subject (native space)
        copyfile(csf_mask,[anat_folder char(org.analysesDate) '/' 'csf_mask.nii' ]);
        csf_mask = p_mvcs_get_files([anat_folder char(org.analysesDate) '/'],'nii','csf_mask'); 
        Vt=spm_vol(csf_mask);
        tmp=spm_read_vols(Vt);
        Vt.pinfo=[0.001;0;0];
        Vt.dt=[8 0];

        if sum(abs(ref_mat(:)-Vt.mat(:)))>0 % if mismatch between affine transform matrices from mask and functional references; will reslice brain mask and over-write image    
            spm_write_vol(Vt,tmp);

            clear matlabbatch;
            matlabbatch{1}.spm.spatial.coreg.write.ref{1} = ref_func;
            matlabbatch{1}.spm.spatial.coreg.write.source{1} = Vt.fname; 
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            spm_jobman('run', matlabbatch);
            clear matlabbatch; 

            Vtx=spm_vol([anat_folder char(org.analysesDate) '/' 'rcsf_mask.nii' ]);
            tmp=spm_read_vols(Vtx);
            Vtx.fname = csf_mask;
            delete([anat_folder char(org.analysesDate) '/csf_mask.nii']);
            delete([anat_folder char(org.analysesDate) '/rcsf_mask.nii']);
            spm_write_vol(Vtx,tmp);
            clear Vtx;   

        end % if sum(abs(ref_mat(:)-Vt.mat(:)))>0 
        tmp=tmp>=prctile(tmp(tmp>0),90); % take those that are greater than 90th percentile 
        mask=tmp(brain_mask_index);
        csf_mask_index=find(mask > 0.5*max(mask(:)));
        clear tmp; clear mask; clear Vt;    
    end
    
    % Gray Matter Mask - used to ensure voxels within ROI are predominantly GM
    if org.preprocessing.gm_check == 1
        disp(['GM mask for subject ' num2str(nSub)]);
        gm_mask = p_mvcs_get_files(anat_folder,'nii', org.identifier_gm); % Find brain mask of subject (native space)
        copyfile(gm_mask,[anat_folder char(org.analysesDate) '/' 'gm_mask.nii' ]);
        gm_mask = p_mvcs_get_files([anat_folder char(org.analysesDate) '/'],'nii','gm_mask'); 
        Vt=spm_vol(gm_mask);
        tmp=spm_read_vols(Vt);
        Vt.pinfo=[0.001;0;0];
        Vt.dt=[8 0];

        if sum(abs(ref_mat(:)-Vt.mat(:)))>0 % if mismatch between affine transform matrices from mask and functional references; will reslice brain mask and over-write image    
            spm_write_vol(Vt,tmp);

            clear matlabbatch;
            matlabbatch{1}.spm.spatial.coreg.write.ref{1} = ref_func;
            matlabbatch{1}.spm.spatial.coreg.write.source{1} = Vt.fname; 
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            spm_jobman('run', matlabbatch);
            clear matlabbatch; 

            Vtx=spm_vol([anat_folder char(org.analysesDate) '/' 'rgm_mask.nii' ]);
            tmp=spm_read_vols(Vtx);
            Vtx.fname = gm_mask;
            delete([anat_folder char(org.analysesDate) '/gm_mask.nii']);
            delete([anat_folder char(org.analysesDate) '/rgm_mask.nii']);
            spm_write_vol(Vtx,tmp);
            clear Vtx;   

        end % if sum(abs(ref_mat(:)-Vt.mat(:)))>0   
        mask=tmp(brain_mask_index);
        gm_mask_index=find(mask > org.preprocessing.gm_check_prob); % voxel must have this probabilty of being GM to be included in analyses
        clear tmp; clear mask; 
        clear Vt; 
    end
    
    % ROI masks
    for nROI = 1:length(org.roi) % loop through ROIs

        disp(['ROI mask ' num2str(nROI) ' for subject ' num2str(nSub)]);
        roi_mask = p_mvcs_get_files(char(anat_folder),'nii', char(org.roi{nROI}));
        xx=org.roi{nROI}; roi.name = xx(2:end); clear xx;
        
        % Get information on ROI
        Vm = spm_vol(roi_mask); 
        mask_roi=spm_read_vols(Vm);
        Vm.fname=[char(anat_folder) char(org.analysesDate) '/cluster_' char(org.roi{nROI}) '.nii'];
        Vm.pinfo=[0.001;0;0];
        Vm.dt=[8 0];
        spm_write_vol(Vm,mask_roi);    
        roi.file = Vm.fname;
        clear mask_roi; clear Vm;
        
        V=spm_vol(roi.file);        
        if sum(abs(ref_mat(:)-V.mat(:)))>0 % if mismatch between affine transform matrices from roi mask and functional reference; will reslice brain mask and over-write image 
            clear matlabbatch;
            matlabbatch{1}.spm.spatial.coreg.write.ref{1} = ref_func;
            matlabbatch{1}.spm.spatial.coreg.write.source{1} = roi.file;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 1;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
            matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
            matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';
            spm_jobman('run', matlabbatch);
            clear matlabbatch;

            delete([char(anat_folder) char(org.analysesDate) '/cluster_w' roi.name '.nii']);
            copyfile([char(anat_folder) char(org.analysesDate) '/rcluster_w' roi.name '.nii'],[char(anat_folder) char(org.analysesDate) '/cluster_w' roi.name '.nii']);
            delete([char(anat_folder) char(org.analysesDate)  '/rcluster_w' roi.name '.nii']);   
            
            ROI = p_mvcs_get_files([char(anat_folder) char(org.analysesDate) '/'],'nii',['cluster_' char(org.roi{nROI})]); 
            
            matlabbatch{1}.spm.util.imcalc.input = {ROI};
            matlabbatch{1}.spm.util.imcalc.output = ['cluster_w' roi.name];
            matlabbatch{1}.spm.util.imcalc.outdir = {[char(anat_folder) char(org.analysesDate)]};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run', matlabbatch);
            clear matlabbatch; clear ROI;
          
        end % if sum(abs(ref_mat(:)-Vm.mat(:)))>0 
        clear V;
        
        ROI = p_mvcs_get_files([char(anat_folder) char(org.analysesDate) '/'],'nii',['cluster_' char(org.roi{nROI})]); 
        % ERODE ROIs
        if org.preprocessing.erodeROI(nROI) == 1
            nii = load_untouch_nii(ROI);
            nii.img = uint8(spm_erode(double(nii.img)));
            save_untouch_nii(nii, ROI); 
            
            ROI = p_mvcs_get_files([char(anat_folder) char(org.analysesDate) '/'],'nii',['cluster_' char(org.roi{nROI})]); 
            
            matlabbatch{1}.spm.util.imcalc.input = {ROI};
            matlabbatch{1}.spm.util.imcalc.output = ['cluster_w' roi.name];
            matlabbatch{1}.spm.util.imcalc.outdir = {[char(anat_folder) char(org.analysesDate)]};
            matlabbatch{1}.spm.util.imcalc.expression = 'i1>0';
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run', matlabbatch);
            clear matlabbatch; 
        end
        
        V=spm_vol(ROI);
        mask_roi=spm_read_vols(V);
        clear ROI; clear V; 
        
        mask_roix=reshape(mask_roi,1,numel(mask_roi)) > 0;
        elem=find(mask_roix(brain_mask_index) > 0);
        if org.preprocessing.gm_check == 1
            roi.elem=intersect(elem,gm_mask_index);
        else
            roi.elem=elem;
        end
        
        if isempty(elem)
            disp(['warning! check the ROI ' roi.file]);
        end % isempty(elem)
        clear elem; clear mask_roix; clear mask_roi; 
        
        % Functional Runs
        for nRun = 1:1:length(org.sessions(nSub,:)) % loop through runs
            if strcmp(org.sessions(nSub,nRun), {'NaN'}) == 0
                disp(['Accessing functional data for subject ' num2str(nSub) '; Run ' num2str(nRun) '; ROI ' num2str(nROI)]);
                clear datat; 
                func_folder = fullfile(org.main_folder{1}, org.subjects{nSub},org.sessions{nSub,nRun}, filesep);
                if not(isfolder([func_folder char(org.analysesDate) '/']))
                    mkdir([func_folder char(org.analysesDate) '/'])
                end

                func_image = p_mvcs_get_files(func_folder,'img', org.identifier_func);
                Vt = spm_vol(func_image); 
                dim(1)=Vt(1).dim(1); dim(2)=Vt(1).dim(2); dim(3)=Vt(1).dim(3); tdim=length(Vt);
                data=spm_read_vols(Vt); % 4D image of functional run data
                data=reshape(data,dim(1)*dim(2)*dim(3),tdim); % reshape such that it is a m x n matrix, where m = voxel and n = volumes 
                datat=data(brain_mask_index,:)';
                clear data;

                % DATA DETRENDING
                datat=detrend(datat); % data detrending is compulsory!          

                % FILTERING
                if org.preprocessing.hp_filter == 1
                    Fs=1/org.TR(nRun);
                    f_high=max(org.preprocessing.hp_value*2/Fs,0.004);
                    [b2,a2] = cheby2(6,20,f_high,'high');
                    disp('high pass filtering...');       
                    datat2=[flipud(datat(1:15,:)) ; datat ; flipud(datat(end-14:end,:))];

                    for zz=1:size(datat,2)
                        datat2(:,zz)=filtfilt(b2,a2,datat2(:,zz));
                    end
                    datat=datat2(16:end-15,:);
                    clear zz; clear datat2; 
                end % if org.preprocessing.hp_filter == 1
                clear Fs; clear f_high; clear b2; clear a2;
                
                % ADJUST FOR DIFFERENT SAMPLING FREQUENCIES (TRs)
                if abs(org.TR_ref-org.TR(nRun))>0
                    disp('adjust for different sampling frequencies...'); 
                    datat=[flipud(datat(1:10,:)) ; datat ; flipud(datat(end-9:end,:))];
                    datat = resample(datat,round(1000*org.TR(nRun)),round(1000*org.TR_ref));  % corrects for different sampling frequencies
                    tdim=ceil(tdim/org.TR_ref*org.TR(nRun));
                    offset=fix((size(datat,1)-tdim)/2);
                    datat=datat(offset+1:offset+tdim,:);
                end % if abs(org.TR_ref-org.TR(nRun))>0
                clear offset; 

                % SCRUBBING
                if org.preprocessing.scrubbing == 1
                    disp('scrubbing...');         

                    motion = p_mvcs_get_files(char(func_folder),'.txt', 'rp_');
                    motion=load(motion,'-ascii');
                    
                    if abs(org.TR_ref-org.TR(nRun))>0
                        motion = resample(motion,round(1000*org.TR(nRun)),round(1000*org.TR_ref));  % corrects motion parameters for different sampling frequencies
                    end % if abs(org.TR_ref-org.TR(nRun))>0
                    
                    val_rot=zeros(tdim,3);
                    val_trans=zeros(tdim,3);
                    for xx=1:tdim
                        x=motion(xx,4); y=motion(xx,5); z=motion(xx,6);
                        val_rot(xx,:)=org.head_size*([1 0 0 ; 0 cos(x) -sin(x); 0 sin(x) cos(x)]*[cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)]*[cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1]-eye(3))*ones(3,1);
                        val_trans(xx,:)=motion(xx,1:3)';
                    end % for xx=1:tdim
                    clear xx; clear x; clear y; clear z;
                    clear fd; 

                    fd=sum(abs(diff(val_trans)),2)+sum(abs(diff(val_rot)),2);  % calculate the framewise displacement
                    fd=[mean(fd) ; fd];  %  accounting for temporal shifts induced by the use of the first derivative                 
                    bad_vols=find(fd >= org.preprocessing.scrubbing_fd_thres); % index of the volumes to be scrubbed
                    if ~isempty(bad_vols)
                        bad_volumes = [];
                        for ii = 1:length(bad_vols)
                            bad_volumes(end+1:end+(org.preprocessing.scrubbing_forwards+org.preprocessing.scrubbing_backwards+1),1) = bad_vols(ii)-org.preprocessing.scrubbing_backwards:1:bad_vols(ii)+org.preprocessing.scrubbing_forwards;
                        end
                        bad_volumes=unique(bad_volumes);   % eliminate redundant bad_volumes
                        bad_volumes=bad_volumes(bad_volumes<=tdim); % eliminate those that are beyond the # of volumes acquired
                        bad_volumes=bad_volumes(bad_volumes>0); % eliminate bad_volume = 0 
                        datat(bad_volumes,:) = NaN;
                        clear ii;
                    else 
                        bad_volumes=bad_vols;
                    end % if ~isempty(bad_vols) 
                    bad_volumes_num = length(bad_volumes); 
                    good_volumes_num = tdim-length(bad_volumes); 
                    bad_volumes_perc = length(bad_volumes)/tdim*100; 
                    clear bad_vols; clear motion; 
                    datat(bad_volumes,:) = [];
                 end % org.preprocessing.scrubbing == 1
                 tdim=length(datat(:,1));

                 warning off;
                 regressors=[];
                 % MOTION REGRESSORS
                 if org.preprocessing.move_regression == 1
                     disp('creating motion regressors...'); 
                     motion = p_mvcs_get_files(char(func_folder),'.txt', 'rp_');
                     motion=load(motion,'-ascii');
                     
                    if abs(org.TR_ref-org.TR(nRun))>0
                        motion = resample(motion,round(1000*org.TR(nRun)),round(1000*org.TR_ref));  % corrects motion parameters for different sampling frequencies
                    end % if abs(org.TR_ref-org.TR(nRun))>0
                    
                    clear motion_shift;
                    motion_shift(size(motion,1),:)=motion(size(motion,1),:);                % changed by DM 04/03/2014
                    motion_shift(1:size(motion,1)-1,:)=motion(2:size(motion,1),:);          % changed by DM 04/03/2014
                    if org.preprocessing.scrubbing == 1
                        motion(bad_volumes,:) = [];
                        motion_shift(bad_volumes,:) = [];
                    end
                    regressors=[motion motion_shift motion.^2 motion_shift.^2 regressors];  % changed by DM 04/03/2014
                    clear motion_shift; clear motion; 
                 end % if org.preprocessing.move_regresssion == 1

                 % WM REGRESSORS
                 if org.preprocessing.wm_regression == 1
                      disp('creating WM regressors...'); 
                      [~, score, ~ , ~ , explained_var.wm] = pca(datat(:,wm_mask_index));
                      wm=score(:,1:org.preprocessing.num_of_pcs);
                      regressors=[wm regressors];
                      clear wm; clear score;
                 end % if org.preprocessing.WM_regresssion == 1

                 % CSF REGRESSORS
                 if org.preprocessing.csf_regression == 1
                      disp('creating CSF regressors...'); 
                       [~, score, ~, ~, explained_var.csf] = pca(datat(:,csf_mask_index));
                       csf=score(:,1:org.preprocessing.num_of_pcs);
                       regressors=[csf regressors];
                       clear csf; clear score; 
                 end % if org.preprocessing.CSF_regresssion == 1
                 
                 % EXTRACT DATA FROM ROI; speeds up removal of regressors
                 datat_roi = datat(:,roi.elem); % data from the ROI
                 roi.voxelposition = voxel_position(:,roi.elem)'; % voxel positions within ROI'
                                 
                 % REMOVE REGRESSORS FROM DATA
                 if ~isempty(regressors)
                       disp('removing regressors from data...');
                       for cont=1:size(datat_roi,2)
                           sig=datat_roi(:,cont);
                           B2=glmfit(regressors,sig,'normal');
                           datat_roi(:,cont)=sig-regressors*B2(2:end);
                       end % for cont=1:length(brain_mask_index)
                       clear cont; clear sig; clear B2;
                 end % if ~isempty(regressors)
                                                 
                 % MATCHING # OF VOLS ACROSS RUNS
                 if org.preprocessing.volume_number_match == 1
                     disp('matching # of vols...');
                     min_volumes = min(precheck.good_volumes_num(nSub,:));
                     if tdim ~= min_volumes
                         datat_roi = datat_roi(round(tdim/2)-floor(min_volumes/2):round(tdim/2)+floor(min_volumes/2),:); % take the middle n volumes, where n is equal to the minimum # of good volumes across all runs for the individual                  
                     end % if tdim ~= min_volumes      
                 end % if org.preprocessing.volume_number_match == 1
                 clear min_volumes;
                 
                % COMPUTE CORRELATIONS AMONG VOXELS WITHIN THE ROI
                clear RHO; clear RHO_zscore;
                RHO = corr(datat_roi);
                RHO_zscore = 0.5 * log((1 + RHO) ./ (1 - RHO));           

                mm=1;
                figure('visible', 'off');
%                 figure(1)
                imagesc(RHO_zscore); axis square; caxis([-mm mm]); colormap(colors); colorbar; daspect([1 1 1]); title(['subj ' num2str(nSub) '; Run ' num2str(nRun) '; ROI ' num2str(nROI)], 'FontSize', 18);
                set(gca,'YTick',[],'YTickLabel',[]);  set(gca,'XTick',[],'XTickLabel',[]);
                print([func_folder char(org.analysesDate) '/pattern_ROI' num2str(nROI) '.tif'],'-dtiff','-r150'); 
                clear mm;
                               
                save([func_folder char(org.analysesDate) '/output_ROI' num2str(nROI) ], 'datat_roi', 'fd', 'regressors', 'RHO_zscore', 'explained_var', 'bad_volumes');                
                
                % CREATE MATRIX WITH ALL RUNS TO COMPUTE SIMILARITY; 
                RHO_zscore_allRuns(:,:,nRun) = RHO_zscore;
                counter = 1;
                for i_col = 1:1:size(RHO_zscore_allRuns,1)-1
                    for i_row = i_col+1:1:size(RHO_zscore_allRuns,2)
                        RHO_zscore_vect_allRuns(counter,nRun) =  RHO_zscore_allRuns(i_row,i_col,nRun);
                        counter = counter + 1;
                    end % for i_row = i_col+1:1:size(RHO_zscore_allRuns,2)
                end % for i_col = 1:1:size(RHO_zscore_allRuns,1)-1
                clear counter; clear i_col; clear i_row;

                scrubbing.bad_volumes_perc_allRuns(nRun) = bad_volumes_perc;
                scrubbing.bad_volumes_num_allRuns(nRun) = bad_volumes_num;
                scrubbing.good_volumes_num_allRuns(nRun) = good_volumes_num;
                
            else % Run is missing (% if strcmp(org.sessions(nSub,nRun), {'NaN'}) == 0)
                
                RHO_zscore_allRuns(1:length(roi.elem),1:length(roi.elem),nRun) = NaN;
                RHO_zscore_vect_allRuns(1:(((length(roi.elem)*length(roi.elem))-length(roi.elem))/2),nRun) =  NaN;
                scrubbing.bad_volumes_perc_allRuns(nRun) = NaN;
                scrubbing.bad_volumes_num_allRuns(nRun) = NaN;
                scrubbing.good_volumes_num_allRuns(nRun) = NaN;
                
            end % if strcmp(org.sessions(nSub,nRun), {'NaN'}) == 0
            
        end % for nRun = 1:1:length(org.sessions(nSub,:))
        Similarity_matrix{1,:,:} = corr(RHO_zscore_vect_allRuns);
        Similarity_matrix_zscore{1,:,:} = 0.5 * log((1 + Similarity_matrix{1,:,:}) ./ (1 - Similarity_matrix{1,:,:})); 
        if not(isfolder([char(org.main_folder) '/' char(org.subjects(nSub)) '/' char(org.analysesDate) '/']))
            mkdir([char(org.main_folder) '/' char(org.subjects(nSub)) '/' char(org.analysesDate) '/'])
        end
        save([char(org.main_folder) '/' char(org.subjects(nSub)) '/' char(org.analysesDate) '/corr_matrices_ROI' num2str(nROI) char(org.analysesDate)], 'Similarity_matrix_zscore', 'scrubbing', 'roi'); 
        clear RHO_zscore_allRuns; clear RHO_zscore_vect_allRuns;
    end % for nROI = 1:length(org.roi) 
end % for nSub =1:length(org.subjects)