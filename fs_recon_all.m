% RKUL_MTLE-HS-fs
% v 1 - dd 18/10/2018
% @ radwan - sunaert
% We will be doing both the T1 and DTI preproc here

% Info:
% This is a simple script that loops over your bias corrected T1s and does recon-all on them.

%% Part 1 Configure my stuff...
clear all;
clc;
dir_main   = '/Users/aradwa0/MR-data/MTLE_HS';
v = ('1 - dd 18/10/2018');

% Other dirs
dir_source_T1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'orig'];
dir_bc_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'bc'];
dir_brains_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'brains'];
dir_output_fs = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'fs_recon-all'];

% need to specify outputs for int. steps of DTI
dir_log    = [dir_output_fs filesep 'LOG'];

% Say hello
str = sprintf('\n%s\n\n','*** Freesurfer recon-all script for MTLE_HS');
disp(str);

% Create folders
% Dir_log should already be created or simply add an mkdir step
% mkdir(dir_log); % create log folder
unix(['echo > ' dir_log filesep 'MTLE_HS_preproc.log']);
mkdir(dir_output_fs);
diary([dir_log filesep 'comm_log_fs_recon-all.txt']);


% Search the source for subjects
s_DTI = (dir([dir_source_DTI filesep '*_1mm.nii.gz']));

% list the T1s found and the DTIs found
DTI_files = (extractfield(s_DTI,'name'))';
T1s_files = (extractfield(s_T1s,'name'))';
unix(['touch ' dir_log filesep 'sanity_check_fail.txt']);
unix(['touch ' dir_log filesep 'sanity_check_pass.txt']);
% quick sanity check to see if all T1s are paired with a DTI
for i = 1:(size(T1s_files,1))
    crs_check = strncmpi(T1s_files(i), DTI_files(i), 6);
    if crs_check == 0
        unix(['echo "*** It seems that ' char(T1s_files(i)) '  has no 1mm T1  damn" >> ' dir_log filesep 'sanity_check_fail.txt' ]);
    else
        unix(['echo "*** It seems that ' char(T1s_files(i)) ' has a 1mm T1   good" >> ' dir_log filesep 'sanity_check_pass.txt' ]);
    end
end

for i = 1:(size(T1s_files,1))
    T1_subj = char(T1s_files(i));
    T1_base_name = (T1_subj(1:end-7));
    T1_subj_bc = ([T1_base_name '_bc.nii.gz']);
    subjid = ([T1_base_name '_fs_rec-all']);
    unix(['recon-all -i ' dir_bc_t1 filesep T1_subj_bc ' -sd ' dir_output_fs filesep  ' -subjid ' subjid '-all -noskullstrip -parallel -openmp 7']);
end

diary off
