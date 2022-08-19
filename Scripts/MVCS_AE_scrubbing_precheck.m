% DIRECTORIES OF INTEREST
org.main_folder = {'/Volumes/DATA/MVCS_ALLO_EGO_20200406/MVCS/DATA_5runs'};
org.analysesDate = {'_analyses_20211022'}; % output data / figures in folder with this Date

% PARAMETERS
org.head_size=50;  % diameter of a standard head (in mm)
org.fd_thres=0.5; % threshold (in mm) for determining bad volumes based on framewise displacement 
org.TR = [2.65; 2.65; 2.65; 2.65; 2.65]; % TR by session (all runs had TR=2.65]
org.TR_ref = 2.65; % Setting the reference TR (resampled such that all data acquired with this TR); recommended to use the longer TR as the default
org.preprocssing.scrubbing = 1; % 0 if step not done; 1 if done
org.preprocssing.scrubbing_backwards = 0; % number of volumes to be removed BEFORE volume detected to be bad
org.preprocssing.scrubbing_forwards = 1;  % number of volumes to be removed AFTER volume detected to be bad
org.preprocssing.volume_number_match = 1; % 0 if step not done; 1 if done; matches the number of volumes across all runs within a participant (after scrubbing and re-sampling for different TR)

files = dir(char(org.main_folder));
counter = 0;
for ii = 1:1:length(files)
    if contains(files(ii).name,'ALLO_EGO')
        counter = counter + 1;
        org.subjects(counter,1) = {files(ii).name};    
    end
end
clear counter; clear ii; clear files;

org.sessions = repmat({'NaN'},length(org.subjects),5);  % 3 is chosen as a max 
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

precheck.bad_volumes_num=NaN(length(org.subjects),5);
precheck.bad_volumes_perc=NaN(length(org.subjects),5);
precheck.good_volumes_num=NaN(length(org.subjects),5);

% MAIN SUBJECT LOOP
for isub = 1:length(org.subjects) % for each subject
    % plot realignment parameters for each subject and session
    for ises = 1:length(org.sessions(isub,:))   
        if strcmp(org.sessions(isub,ises), {'NaN'}) == 0
            func_folder = fullfile(org.main_folder{1}, org.subjects{isub},org.sessions{isub,ises}, filesep);
            motion = p_mvcs_get_files(char(func_folder),'.txt', 'rp_');
            motion=load(motion,'-ascii');
            
            if abs(org.TR_ref-org.TR(ises))>0
                motion = resample(motion,round(1000*org.TR(ises)),round(1000*org.TR_ref));  % corrects motion parameters for different sampling frequencies
            end % if abs(org.TR_ref-org.TR(nRun))>0
            tdim=length(motion(:,1));
            
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
            bad_vols=find(fd >= org.fd_thres); % index of the volumes to be scrubbed
            if ~isempty(bad_vols)
                bad_volumes = [];
                for ii = 1:length(bad_vols)
                    bad_volumes(end+1:end+(org.preprocssing.scrubbing_forwards+org.preprocssing.scrubbing_backwards+1),1) = bad_vols(ii)-org.preprocssing.scrubbing_backwards:1:bad_vols(ii)+org.preprocssing.scrubbing_forwards;
                end
                bad_volumes=unique(bad_volumes);   % eliminate redundant bad_volumes
                bad_volumes=bad_volumes(bad_volumes<=tdim); % eliminate those that are beyond the # of volumes acquired
                bad_volumes=bad_volumes(bad_volumes>0); % eliminate bad_volume = 0 
                clear ii;
            else 
                bad_volumes=bad_vols;
            end % if ~isempty(bad_vols) 
            precheck.bad_volumes_num(isub,ises) = length(bad_volumes); 
            precheck.good_volumes_num(isub,ises) = tdim-length(bad_volumes); 
            precheck.bad_volumes_perc(isub,ises) = length(bad_volumes)/tdim*100; 
            clear bad_vols; clear motion; 
        end

    end    
end
save([char(org.main_folder) '/_scrubbing_precheck/'  'scrubbing_precheck' char(org.analysesDate)], 'org', 'precheck'); 
