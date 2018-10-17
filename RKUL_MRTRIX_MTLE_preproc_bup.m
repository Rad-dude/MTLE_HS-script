% RKUL_MTLE-HS-preproc
% v 1.7 - dd 15/10/2018
% @ radwan - sunaert
% We will be doing both the T1 and DTI preproc here

% Info:
% Here we're starting the prepped and anonimyzed nii T1s and DTIS,
% preprocessing with ANTs for T1s and with ANTS, FSL and MRTRIX for the
% DTIs
% so far (13/10) the script finishes off by producing prob_tensor whole
% brain tractograms with 5 million streamlines (using act and seed_gmwmi),
% and the AAL and BN246 labels in native DTI space.
% We rely on nonlinear registration based EDC for the DTIs and we're using
% antsBET and antsAtropos directly on the B0s (more reliable than dwi2mask
% and 5ttgen)

%% Part 1 Configure my stuff...
clear all;
clc;
dir_main   = '/Users/aradwa0/MR-data/MTLE_HS';
v = ('1.7 - dd 15/10/2018');

% Other dirs
% dir_source_DTI = [dir_main filesep 'NII' filesep 'DTI' filesep 'orig'];
% dir_source_T1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'orig'];
% dir_bc_t1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'bc'];
% dir_brains_t1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'brains'];
% dir_tpms_t1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'tpms'];
% dir_labels_t1 = [dir_main filesep 'NII' filesep 'T1s' filesep 'labels'];
% dir_output_DTI = [dir_main filesep 'NII' filesep 'DTI' filesep 'MIF'];

dir_source_DTI = [dir_main filesep 'NII_all' filesep 'DTI' filesep 'orig'];
dir_source_T1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'orig'];
dir_bc_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'bc'];
dir_brains_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'brains'];
dir_tpms_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'tpms'];
dir_labels_t1 = [dir_main filesep 'NII_all' filesep 'T1s' filesep 'labels'];
dir_output_DTI = [dir_main filesep 'NII_all' filesep 'DTI' filesep 'MIF'];

dir_DTI_pping = [dir_output_DTI filesep 'preprocessing']; % dir for intermediate DTI outputs
dir_DTI_ready = [dir_output_DTI filesep 'DTI_ready']; % dir for prepped DTI outputs

% need to specify outputs for int. steps of DTI
% dir_T1_reconall = [dir_main filesep 'fs_rec_all'];
dir_log    = [dir_main filesep 'LOG'];
dir_templates = [dir_main filesep 'templates'];

% Say hello
str = sprintf('\n%s\n\n','*** DTI preproc script for MTLE_HS');
disp(str);

% Create folders
% Dir_log should already be created or simply add an mkdir step
% mkdir(dir_log); % create log folder
unix(['echo > ' dir_log filesep 'MTLE_HS_preproc.log']);
mkdir(dir_output_DTI);
mkdir(dir_bc_t1);
mkdir(dir_brains_t1);
mkdir(dir_tpms_t1);
mkdir(dir_labels_t1);
mkdir(dir_DTI_pping);
mkdir(dir_DTI_ready);
% mkdir(dir_T1_reconall);
diary([dir_log filesep 'command_log_dwi_preproc.txt']);

% define templates
T1_template = [dir_templates filesep 'PD25-T1MPRAGE-template-1mm.nii.gz'];
T1_template_brain = [dir_templates filesep 'PD25-T1MPRAGE-template-1mm_brain.nii.gz'];
T1_temp_brain_mask = [dir_templates filesep 'PD25-atlas-mask-1mm.nii.gz'];
ttpm_gm_cx = [dir_templates filesep 'PD25_priors1.nii.gz'];
ttpm_gm_bg = [dir_templates filesep 'PD25_priors2.nii.gz'];
ttpm_wm = [dir_templates filesep 'PD25_priors3.nii.gz'];
ttpm_csf = [dir_templates filesep 'PD25_priors4.nii.gz'];
aal_template_brain = [dir_templates filesep 'MNI152_T1_1mm_brain.nii.gz'];
aal_labels = [dir_templates filesep 'aal_atlas.nii'];
txt_aal_labels = [dir_templates filesep 'aal_atlas.txt'];
BN_246_3D_labels = [dir_templates filesep 'BN_Atlas_246_1mm.nii.gz'];
% BN_246_template_brain = [dir_templates filesep 'PD25-T1MPRAGE-template-1mm_brain.nii.gz'];
txt_BN_246_3D_labels = [dir_templates filesep 'BN_Atlas_246_LUT.txt'];
T2_template = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a.nii'];
T2_temp_brain_mask = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a_mask.nii'];
T2_temp_brain = [dir_templates filesep 'mni_icbm152_t2_tal_nlin_asym_09a_brain.nii.gz'];


% Search the source for subjects
s_DTI = (dir([dir_source_DTI filesep '*.nii.gz']));
s_T1s = (dir([dir_source_T1 filesep '*_orig.nii.gz']));

% list the T1s found and the DTIs found
DTI_files = (extractfield(s_DTI,'name'))';
T1s_files = (extractfield(s_T1s,'name'))';
unix(['touch ' dir_log filesep 'sanity_check_fail.txt']);
unix(['touch ' dir_log filesep 'sanity_check_pass.txt']);
% quick sanity check to see if all T1s are paired with a DTI
for i = 1:(size(T1s_files,1))
    crs_check = strncmpi(T1s_files(i), DTI_files(i), 6);
    if crs_check == 0
        unix(['echo "*** It seems that ' char(T1s_files(i)) ' has no DWI data" >> ' dir_log filesep 'sanity_check_fail.txt' ]);
    else
        unix(['echo "*** It seems that ' char(T1s_files(i)) ' has DWI data" >> ' dir_log filesep 'sanity_check_pass.txt' ]);
    end
end

%% Part 2 - T1s preproc with ANTs

for i = 1:(size(T1s_files,1))
    T1_subj = char(T1s_files(i));
    T1_base_name = (T1_subj(1:end-7));
    % something else that should be moved to a struct rather than defined with
    % every loop.
    T1_subj_1mm = ([T1_base_name '_1mm.nii.gz']);
    T1_subj_bc = ([T1_base_name '_bc.nii.gz']);
    T1_subj_bf = ([T1_base_name '_bf.nii.gz']);
    T1_subj_brain = ([T1_base_name '_']);
    Temp_brain2subj_brain = (['template_brain_2_' (T1_subj_brain(1:end-13))]);
    aal_in_subj = [dir_labels_t1 filesep 'aal_in' (T1_base_name(1:end-12)) '.nii.gz'];
    BN246_in_subj = [dir_labels_t1 filesep 'BN246_in' (T1_base_name(1:end-12)) '.nii.gz'];
    ttpm_gm_cx_in_subj = (['PD25-GM_cx_in_' (T1_base_name(1:end-12)) '_brain_1.nii.gz']);
    ttpm_gm_bg_in_subj = (['PD25-GM_bg_in_' (T1_base_name(1:end-12)) '_brain_2.nii.gz']);
    ttpm_wm_in_subj = (['PD25-WM_in_' (T1_base_name(1:end-12)) '_brain_3.nii.gz']);
    ttpm_csf_in_subj = (['PD25-CSF_in_' (T1_base_name(1:end-12)) '_brain_4.nii.gz']);
    %     First we reslice the T1s to 1 mm isotropic
    unix(['source ~/.bash_profile ; mrresize -force -nthreads 7 -voxel 1 ' dir_source_T1 filesep T1_subj ' ' dir_source_T1 filesep T1_subj_1mm]);
    unix(['N4BiasFieldCorrection -d 3 -i ' dir_source_T1 filesep T1_subj_1mm ...
        ' -o [' dir_bc_t1 filesep T1_subj_bc ',' dir_bc_t1 filesep T1_subj_bf ']']);
    % alternative is to use ANTS yo!
    % BET
    unix(['antsBrainExtraction.sh -d 3 -a ' dir_bc_t1 filesep T1_subj_bc ' -e ' T1_template ' -m ' T1_temp_brain_mask ...
        ' -o ' dir_brains_t1 filesep T1_subj_brain ' -s nii.gz -u 1']);
    % Grab template and priors to subj space
    unix(['antsRegistrationSyNQuick.sh -d 3 -m ' T1_template_brain ...
        ' -f ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz -o ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain ' -x ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -n 7 -j 1 -t s']);
    unix(['WarpImageMultiTransform 3 ' ttpm_csf ' ' dir_tpms_t1 filesep ttpm_csf_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' ttpm_wm ' ' dir_tpms_t1 filesep ttpm_wm_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' ttpm_gm_cx ' ' dir_tpms_t1 filesep ttpm_gm_cx_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' ttpm_gm_bg ' ' dir_tpms_t1 filesep ttpm_gm_bg_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
        dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    % T1_brain segmentation into Cortex, BG, WM,& CSF
    unix(['antsAtroposN4.sh -d 3 -a ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz -x ' ...
        dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -c 4 -p ' ...
        dir_tpms_t1 filesep (T1_base_name(1:end-12)) '%d.nii.gz -r 0.5 -w 0.05 -o ' dir_tpms_t1 filesep T1_base_name 'atropos_' ]);
    % Warping Atlas templates and labels
    unix(['antsRegistrationSyNQuick.sh -d 3 -m ' aal_template_brain ...
        ' -f ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz -o ' ...
        dir_labels_t1 filesep 'MNI1mbrain_2_' (T1_base_name(1:end-12)) 'brain  -x ' ...
        dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -n 7 -j 1 -t s']);
    unix(['antsRegistrationSyNQuick.sh -d 3 -m ' T1_template_brain ... % there was a typo here, now it should be okay 5/10/2018 R
        ' -f ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz -o ' ...
        dir_labels_t1 filesep 'PD25brain_2_' (T1_base_name(1:end-12)) 'brain  -x ' ...
        dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -n 7 -j 1 -t s']);
    unix(['WarpImageMultiTransform 3 ' aal_labels ' ' aal_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_labels_t1 filesep 'MNI1mbrain_2_' (T1_base_name(1:end-12)) 'brain1Warp.nii.gz ' dir_labels_t1 filesep 'MNI1mbrain_2_' (T1_base_name(1:end-12)) 'brain0GenericAffine.mat --use-NN']); % still using NN as ML is rubbish
    unix(['WarpImageMultiTransform 3 ' BN_246_3D_labels ' ' BN246_in_subj ...
        ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
        dir_labels_t1 filesep 'PD25brain_2_' (T1_base_name(1:end-12)) 'brain1Warp.nii.gz ' dir_labels_t1 filesep 'PD25brain_2_' (T1_base_name(1:end-12)) 'brain0GenericAffine.mat --use-NN']);
    
    % will be using the BN_246 labels atlas ? or just AAL ?
    % BN_246 is better but more strict, higher chance of missing some
    % streamlines, AAL is less cool, also larger voxels generally,
    % requiring warping, and lacks a template to use for this purpose!
    % BOTH man!! come on keep up.
    % 04/10/2018
    % For some reason directlabels.sh keeps failing!!
    %     unix(['directlabels.sh -d 3 -s ' dir_tpms_t1 filesep 'atropos_Segmentation.nii.gz -l ' aal_in_subj ...
    %         ' -o aal_in' T1_base_name 'thickness.nii.gz -g ' dir_tpms_t1 filesep ttpm_gm_in_subj (T1_base_name(1:end-12)) '.nii.gz  -w ' ...
    %         dir_tpms_t1 filesep ttpm_wm_in_subj (T1_base_name(1:end-12)) '.nii.gz -c=2 -j 4']); % the input for this argument is actually the segmentation 3D image resulting from Atropos
    %     unix(['directlabels.sh -d 3 -s ' dir_tpms_t1 filesep 'atropos_Segmentation.nii.gz -l ' BN246_in_subj ...
    %         ' -o BN246_in' T1_base_name 'thickness.nii.gz -g ' dir_tpms_t1 filesep ttpm_gm_in_subj (T1_base_name(1:end-12)) '.nii.gz  -w ' dir_tpms_t1 filesep ttpm_wm_in_subj ...
    %         '- c=2 -j 4']);
    
    %     unix(['recon-all -i ' dir_source_DTI filesep T1_subj_bc ' -s ' ...
    %         dir_T1_reconall filesep T1_base_name 'recon_all_output' ' -all']);
    % need to add step for mri_convert brain.mgz to brain.nii
    % fslmaths -thr 10 -bin brainmask.nii
    
    % end
    
    
    
    
    %% Part 3 - DTI_preproc
    
    % for i = 1:(size(DTI_files,1));
    % consider redefining all these loop variables in a struct
    % and just use the index to flip through them.
    
    % Initial variables
    DTI_sub_base = char(DTI_files(i)); % setting the base name
    DTI_sub_base = (DTI_sub_base(1:end-7));
    DTI_bval = [DTI_sub_base '.bval']; % diff weighting input
    DTI_bvec = [DTI_sub_base '.bvec'];
    DTI_bmatrix = [DTI_sub_base '_bmatrix.txt']; % diff weighting int.
    DTI_sub_nii = [DTI_sub_base '.nii.gz'];
    DTI_sub_mif = [DTI_sub_base '.mif'];
    
    % Resized DTI
    upsamp_DTI_sub_mif = [DTI_sub_base '_1mm.mif'];
    
    % Step 1 - denoising
    DTI_1 = [DTI_sub_base '_dn.mif'];
    
    % Step 2 - degibbs
    DTI_2 = [DTI_sub_base '_dn_dg.mif'];
    
    % Step 3 - dwipreproc
    DTI_3 = [DTI_sub_base '_pp.mif'];
    
    % Step 4 - Bias field correction
    DTI_4 = [DTI_sub_base '_pp_bc.mif']; % bias correct
    DTI_4_nii = [DTI_sub_base '_pp_bc.nii']; % semi_preped dti nii
    DTI_4_mif = [DTI_sub_base '_pp_bc.mif']; % semi_preped dti nii
    DTI_4_DT = [DTI_sub_base '_pp_bc_DT.mif']; % diffusion tensor mif
    DTI_4_DT_nii = [DTI_sub_base '_pp_bc_DT.nii'];
    DTI_4_b0 = [DTI_sub_base '_pp_bc_b0.nii'];
    
    % Step 5 - prep for BET of DWI
    b0_brain_basename = [DTI_sub_base '_b0_'];
    DTI_mask_1 = [DTI_sub_base '_pp_bc_mask1.nii'];
    DTI_BET_mask = [b0_brain_basename 'BrainExtractionMask.nii.gz'];
    DTI_FA_4 = [DTI_sub_base '_4_FA_orig.nii'];
    DTI_FA_5 = [DTI_sub_base '_5_FA_BET.nii'];
    DTI_ADC_5 = [DTI_sub_base '_5_ADC_BET.nii'];
    
    % Step 6 - prep for SyN based EDC
    DTI_T1_6 = [DTI_sub_base '_T1_2_FA_affdeformed.nii.gz'];
    DTI_T1_6_mask = [DTI_sub_base '_T1_2_FA_affdeformed_mask.nii.gz'];
    T1_sub_in_FA_lin = [DTI_sub_base '_T1wskull_2_FA_affdeformed.nii.gz'];
    
    % Step 7 - SyN based EDC
    DTI_7 = [DTI_sub_base '_7_edc1.nii.gz'];
    DTI_7_mif = [DTI_sub_base '_7_edc1.mif'];
    DTI_7_DT_nii = [DTI_sub_base '_7_DT_edc1.nii.gz'];
    DTI_7_DT_mif = [DTI_sub_base '_7_DT_edc1.mif'];
    DTI_b0_edc_1_wskull = [DTI_sub_base '_edc1_b0_wskull.nii.gz'];
    
    % Step 8 - BET after EDC
    DTI_BET_mask_edc = [DTI_sub_base '_edc1_BrainExtractionMask.nii.gz']; % still needs editing
    DTI_FA_13 = [DTI_sub_base '_13_DT_postDWISyN.nii.gz'];
    DTI_b0_edc_1 = [DTI_sub_base '_edc1_b0.nii.gz'];
    DTI_b0_edc_1_brain_mask = [b0_brain_basename 'edc1_BrainExtractionMask.nii.gz'];
    DTI_b0_edc_1_brain = [b0_brain_basename 'edc1_BrainExtractionBrain.nii.gz'];
    DTI_FA_edc_1 = [DTI_sub_base '_edc1_FA.nii.gz'];
    DTI_ADC_edc_1 = [DTI_sub_base '_edc1_ADC.nii.gz'];
    
    % Step 9 - Dual channel segmentation of T1 and B0
    t2tpm_gm_cx = [dir_templates filesep 'mni_T2_template_priors1.nii.gz'];
    t2tpm_gm_bg = [dir_templates filesep 'mni_T2_template_priors2.nii.gz'];
    t2tpm_wm = [dir_templates filesep 'mni_T2_template_priors3.nii.gz'];
    t2tpm_csf = [dir_templates filesep 'mni_T2_template_priors4.nii.gz'];
    t2tpm_gm_cx_in_subj = (['T2_temp_in_' DTI_sub_base '_b0_p1.nii.gz']);
    t2tpm_gm_bg_in_subj = (['T2_temp_in_' DTI_sub_base '_b0_p2.nii.gz']);
    t2tpm_wm_in_subj = (['T2_temp_in_' DTI_sub_base '_b0_p3.nii.gz']);
    t2tpm_csf_in_subj = (['T2_temp_in_' DTI_sub_base '_b0_p4.nii.gz']);
    b0_gm_cx = ([DTI_sub_base '_edc1_b0_atropos_SegmentationPosteriors1.nii.gz']);
    b0_gm_bg = ([DTI_sub_base '_edc1_b0_atropos_SegmentationPosteriors2.nii.gz']);
    b0_wm = ([DTI_sub_base '_edc1_b0_atropos_SegmentationPosteriors3.nii.gz']);
    b0_csf = ([DTI_sub_base '_edc1_b0_atropos_SegmentationPosteriors4.nii.gz']);
    b0_empty_image = ([DTI_sub_base '_edc1_b0_p5.nii.gz']); % only to conform to 5tt
    b0_pseudo_5ttgen = ([DTI_sub_base '_edc1_b0_pseudo_5tt.nii.gz']);
    failed_5tt = ([DTI_sub_base '_edc1_b0_p5tt_fail_masks.nii.gz']);
    b0_p5tt_gmwmi = ([DTI_sub_base '_edc1_b0_p5tt_gmwmi.nii.gz']);
    subj_wbrain_prob_tck = ([DTI_sub_base '_edc1_act_gwmwi_probT.tck']);
    subj_wbrain_msmt_csd_tck = ([DTI_sub_base '_edc1_act_gwmwi_msmt_csd.tck']);
    sift_probT_tck = ([DTI_sub_base '_edc1_act_gwmwi_probT_sifted.tck']);
    sift_msmt_csd_tck = ([DTI_sub_base '_edc1_act_gwmwi_msmt_csd_sifted.tck']);
    aal_in_subj_DTI = (['aal_in' DTI_sub_base '.nii.gz']);
    BN246_in_subj_DTI = (['BN246-in' DTI_sub_base '.nii.gz']);
    
    %     % to segment with the T1 as the primary - note to self: this seems
    %     % redundant now, as both T1 and B0 segmentations have the same problem.
    %     % This problem being that they both consider the superior aspect of the
    %     % external capsule as well as the extreme capsule to be mostly GM, how
    %     % will ACT and seedgmwmi deal with this ? Do we lose all sensitivity to
    %     % fibers in this particular part of the brain ?
    %     ttpm_gm_cx_in_dti = (['PD25-GM_cx_in_' DTI_sub_base '_brain_1.nii.gz']);
    %     ttpm_gm_bg_in_dti = (['PD25-GM_bg_in_' DTI_sub_base '_brain_2.nii.gz']);
    %     ttpm_wm_in_dti = (['PD25-WM_in_' DTI_sub_base '_brain_3.nii.gz']);
    %     ttpm_csf_in_dti = (['PD25-CSF_in_' DTI_sub_base '_brain_4.nii.gz']);
    
    %     % Get dims of DTI
    %     [a,b] = unix(['source ~/.bash_profile ; mrinfo -spacing ' dir_DTI_pping filesep DTI_FA_4]);
    %     B = textscan(b,'%s','Delimiter',' ')';
    %     e = char(B{1}{1});
    %     e = e(1:end-13);
    %     f = char(B{1}{2});
    %     f = f(1:end-13);
    %     g = char(B{1}{3});
    %     g = g(1:end-13);
    %     orig_DTI_vox = [e ',' f  ',' g];
    
    % Preparing the DTI data
    unix(['source ~/.bash_profile ; mrconvert -nthreads 7 -force -fslgrad ' dir_source_DTI filesep DTI_bvec ' ' dir_source_DTI filesep DTI_bval ...
        ' ' dir_source_DTI filesep DTI_sub_nii ' ' dir_DTI_pping filesep DTI_sub_mif ' -export_grad_mrtrix ' dir_DTI_pping filesep DTI_bmatrix]);
    % Denoise and degibbs
    unix(['source ~/.bash_profile ; dwidenoise -force -nthreads 7 ' dir_DTI_pping filesep DTI_sub_mif ' ' dir_DTI_pping filesep DTI_1]); % step 1 denoise
    unix(['source ~/.bash_profile ; mrdegibbs -force -nthreads 7 ' dir_DTI_pping filesep DTI_1 ' ' dir_DTI_pping filesep DTI_2]);
    % step 2 - degibbs to correct for gibb's ringing, 2) should add axes 1,2 is for sagittal
    % Upsample to 1 mm isotropic
    unix(['source ~/.bash_profile ; mrresize -force -nthreads 7 -voxel 1 ' dir_DTI_pping filesep DTI_sub_mif ' ' dir_DTI_pping filesep upsamp_DTI_sub_mif]);
    % continue preproc
    unix(['source ~/.bash_profile ; time dwipreproc -force -nthreads 7 -rpe_none -pe_dir AP -eddy_options " --repol --residuals" ' ...
        dir_DTI_pping filesep upsamp_DTI_sub_mif ' ' dir_DTI_pping filesep DTI_3]); % added time ahead of dwipreproc for timing this step
    unix(['source ~/.bash_profile ; dwibiascorrect -force -ants -nthreads 6 ' dir_DTI_pping filesep DTI_3 ' ' dir_DTI_pping filesep DTI_4]);
    unix(['source ~/.bash_profile ; mrconvert -force -nthreads 7 -export_grad_mrtrix ' dir_DTI_pping filesep DTI_bmatrix ' ' dir_DTI_pping filesep DTI_4 ...
        ' ' dir_DTI_pping filesep DTI_4_nii]);
    unix(['source ~/.bash_profile ; dwi2mask -force -nthreads 7 ' dir_DTI_pping filesep DTI_4 ' ' dir_DTI_pping filesep DTI_mask_1]);% step 3
    % No premasking yet
    unix(['source ~/.bash_profile ; dwi2tensor -force -nthreads 7 ' dir_DTI_pping filesep DTI_4 ' ' dir_DTI_pping filesep DTI_4_DT]);
    unix(['source ~/.bash_profile ; dwiextract -force -bzero -nthreads 7 ' dir_DTI_pping filesep DTI_4 ' ' dir_DTI_pping filesep  DTI_4_b0]);
    % premasking initial tensor and FA
    unix(['source ~/.bash_profile ; dwi2tensor -force -nthreads 7 -mask ' dir_DTI_pping filesep DTI_mask_1 ' ' dir_DTI_pping filesep DTI_4 ' ' dir_DTI_pping filesep DTI_4_DT_nii]);% also will need a .nii tensor image
    
    % % % Use this if you want to do a dual channel BET
    % %         unix(['source ~/.bash_profile ; tensor2metric -force -nthreads 7 -mask ' dir_DTI_pping filesep DTI_mask_1 ' -fa ' dir_DTI_pping filesep DTI_FA_4 ...
    % %             ' ' dir_DTI_pping filesep DTI_4_DT]);
    % %     unix(['antsRegistrationSyNQuick.sh -d 3 -m ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz '  ...
    % %         ' -f ' dir_DTI_pping filesep DTI_FA_4 ' -o ' dir_DTI_pping filesep DTI_sub_base '_T1_2_FA_syn -x ' ...
    % %         dir_DTI_pping filesep DTI_mask_1 ' ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -n 7 -j 1 -t s']);
    % % %     need to apply warp to the non-bet T1 WIs using the masks, then do a
    % % %     dual channel brain extraction.. same can be then done also for the
    % % %     Atropos segmentation of the B0
    % %     unix(['WarpImageMultiTransform 3 ' dir_brains_t1 filesep T1_subj ' ' dir_DTI_pping filesep 'T1_sub_in_FA.nii.gz' ...
    % %         ' -R ' dir_DTI_pping filesep DTI_FA_5 ' ' ...
    % %         dir_DTI_pping filesep DTI_sub_base 'T1_2_FA_syn11Warp.nii.gz ' ...
    % %         dir_DTI_pping filesep DTI_sub_base 'T1_2_FA_syn10GenericAffine.mat --use-BSpline']);
    
    % extracting B0 brain
    unix(['antsBrainExtraction.sh -d 3 -a ' dir_DTI_pping filesep DTI_4_b0 ' -e ' T2_template ' -m ' T2_temp_brain_mask ...
        ' -o ' dir_DTI_pping filesep b0_brain_basename ' -s nii.gz -u 1']);
    % regenerating FA with new B0 brain mask
    unix(['source ~/.bash_profile ; tensor2metric -force -nthreads 7 -mask ' dir_DTI_pping filesep DTI_BET_mask ' -fa ' dir_DTI_pping filesep DTI_FA_5 ' -adc ' dir_DTI_pping filesep DTI_ADC_5 ...
        ' ' dir_DTI_pping filesep DTI_4_DT]);
    %     % now brain extraction works well for the B0... the problem was the
    %     % voxels with high values out in the air resulting from dwi2tensor.
    %
    %     % with reg mask
    % %     unix(['antsBrainExtraction.sh -d 3 -a ' dir_DTI_pping filesep DTI_b0 ' -e ' T2_template ' -m ' T2_temp_brain_mask ...
    % %         ' -o ' dir_DTI_pping filesep b0_brain_basename ' -f ' dir_DTI_pping filesep DTI_mask_1 ' -s nii.gz -u 1 -k 1']);
    %
    % %     unix(['antsBrainExtraction.sh -d 3 -a ' dir_DTI_pping filesep DTI_b0 ' -e ' T1_template ' -m ' T1_temp_brain_mask ...
    % %         ' -o ' dir_DTI_pping filesep b0_brain_basename ' -s nii.gz -u 1 -k 1']);
    % %     unix(['antsBrainExtraction.sh -d 3 -a ' dir_DTI_pping filesep DTI_b0 ' -e ' T2_template ' -m ' T2_temp_brain_mask ...
    % %         ' -o ' dir_DTI_pping filesep 'trial_b0_bet_noupsamp -s nii.gz -u 1 -k 1']);
    %     % need to downsample after BET
    % %     unix(['fslmaths ' dir_DTI_pping filesep upsamp_FA_5 ' -mas ' dir_DTI_pping filesep upsamp_DTI_mask ' ' dir_DTI_pping filesep DTI_FA_6 ]);
    % %     unix(['source ~/.bash_profile ; mrresize -voxel ' orig_DTI_vox ' ' dir_DTI_pping filesep DTI_FA_6 ' ' dir_DTI_pping filesep DTI_FA_7]);
    % %     unix(['source ~/.bash_profile ; mrresize -voxel ' orig_DTI_vox ' ' dir_DTI_pping filesep DTI_FA_5 ' ' dir_DTI_pping filesep DTI_FA_7_mask]);
    % %     unix(['fslmaths ' dir_DTI_pping filesep DTI_FA_7_mask ' -bin ' dir_DTI_pping filesep DTI_FA_7_mask]);
    %
    %% Here is the most interesting part (SyN based EDC)
    
    % Linear warp T1 2 FA
    unix(['antsaffine.sh 3 ' dir_DTI_pping filesep DTI_FA_5 ' ' ...
        dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' dir_DTI_pping filesep DTI_sub_base '_T1_2_FA_aff PURELY-RIGID']);
    % Bring your orig_T1_wskull to FA
    unix(['WarpImageMultiTransform 3 ' dir_source_T1 filesep T1_subj_1mm ' ' dir_DTI_pping filesep T1_sub_in_FA_lin ...
        ' -R ' dir_DTI_pping filesep DTI_FA_5 ' ' ...
        dir_DTI_pping filesep DTI_sub_base '_T1_2_FA_affAffine.txt --use-BSpline']);
    % maybe we need to make a new antsRegistrationSyN.sh version with -g,
    % --restrict-deformation PxQxR to restrict deformation to AP axis only.
    unix(['fslmaths ' dir_DTI_pping filesep DTI_T1_6 ' -bin ' dir_DTI_pping filesep DTI_T1_6_mask ]);% binarize T1_in_FA
    unix(['antsRegistrationSyNQuick.sh -d 3 -f ' dir_DTI_pping filesep DTI_T1_6 ' -m ' dir_DTI_pping filesep DTI_FA_5 ' -o ' ...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN -x ' ...
        dir_DTI_pping filesep DTI_T1_6_mask ' ' dir_DTI_pping filesep DTI_BET_mask ' -n 7 -j 1 -t s']);% proceeding with 3 stage Rig,Aff,SyN
    % Bring B0_wskull to T1_in_FA
    unix(['WarpImageMultiTransform 3 ' dir_DTI_pping filesep DTI_4_b0 ' ' dir_DTI_pping filesep DTI_b0_edc_1_wskull ...
        ' -R ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyNWarped.nii.gz '...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN1Warp.nii.gz ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN0GenericAffine.mat --use-BSpline']);
    % Warping the orig_DTI data  % --- testing
    unix(['WarpTimeSeriesImageMultiTransform 4 ' dir_DTI_pping filesep DTI_4_nii ' ' dir_DTI_pping filesep DTI_7 ' -R ' ...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyNWarped.nii.gz '...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN1Warp.nii.gz ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN0GenericAffine.mat --use-BSpline']);
    unix(['source ~/.bash_profile ; mrconvert -force -nthreads 7 -grad ' dir_DTI_pping filesep DTI_bmatrix ' ' dir_DTI_pping filesep DTI_7 ...
        ' ' dir_DTI_pping filesep DTI_7_mif]);
    unix(['source ~/.bash_profile ; dwiextract -force -bzero -nthreads 7 ' dir_DTI_pping filesep DTI_7_mif ' ' dir_DTI_pping filesep  DTI_b0_edc_1]);
    unix(['antsBrainExtraction.sh -d 3 -a ' dir_DTI_pping filesep DTI_b0_edc_1 ' -a ' dir_DTI_pping filesep T1_sub_in_FA_lin ' -e ' T2_template ' -m ' T2_temp_brain_mask ...
        ' -o ' dir_DTI_pping filesep b0_brain_basename 'edc1_ -s nii.gz -u 1']);
    unix(['source ~/.bash_profile ; dwi2tensor -force -nthreads 7 -mask ' dir_DTI_pping filesep DTI_b0_edc_1_brain_mask ...
        ' ' dir_DTI_pping filesep DTI_7_mif ' ' dir_DTI_pping filesep DTI_7_DT_nii]); % consider doing DTIFIT and correcting the tensors...
    unix(['source ~/.bash_profile ; dwi2tensor -force -nthreads 7  -mask ' dir_DTI_pping filesep DTI_b0_edc_1_brain_mask ' ' ...
        dir_DTI_pping filesep DTI_7_mif ' ' dir_DTI_pping filesep DTI_7_DT_mif]);
    unix(['source ~/.bash_profile ; tensor2metric -force -nthreads 7  -fa ' dir_DTI_pping filesep DTI_FA_edc_1 ' -adc ' dir_DTI_pping filesep DTI_ADC_edc_1 ...
        ' ' dir_DTI_pping filesep DTI_7_DT_mif]);
    % Bring template T2 to native_edc_B0 for direct B0 segmentation
    unix(['antsRegistrationSyNQuick.sh -d 3 -m ' T2_temp_brain ' -f ' dir_DTI_pping filesep DTI_b0_edc_1_brain ...
        ' -o ' dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base ' -x ' ...
        dir_DTI_pping filesep DTI_b0_edc_1_brain_mask ' ' T2_temp_brain_mask ' -n 7 -j 1 -t s']);
    unix(['WarpImageMultiTransform 3 ' t2tpm_gm_cx ' ' dir_DTI_pping filesep t2tpm_gm_cx_in_subj ...
        ' -R ' dir_DTI_pping filesep DTI_b0_edc_1 ' ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '1Warp.nii.gz ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' t2tpm_gm_bg ' ' dir_DTI_pping filesep t2tpm_gm_bg_in_subj ...
        ' -R ' dir_DTI_pping filesep DTI_b0_edc_1 ' ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '1Warp.nii.gz ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' t2tpm_wm ' ' dir_DTI_pping filesep t2tpm_wm_in_subj ...
        ' -R ' dir_DTI_pping filesep DTI_b0_edc_1 ' ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '1Warp.nii.gz ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '0GenericAffine.mat --use-BSpline']);
    unix(['WarpImageMultiTransform 3 ' t2tpm_csf ' ' dir_DTI_pping filesep t2tpm_csf_in_subj ...
        ' -R ' dir_DTI_pping filesep DTI_b0_edc_1 ' ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '1Warp.nii.gz ' ...
        dir_DTI_pping filesep 't2_template_2_B0_SyN_' DTI_sub_base '0GenericAffine.mat --use-BSpline']);
    % Dual channel Atropos segmentation of B0 & T1
    unix(['antsAtroposN4.sh -d 3 -a ' dir_DTI_pping filesep DTI_b0_edc_1_brain ' -a ' dir_DTI_pping filesep DTI_T1_6 ' -x ' ...
        dir_DTI_pping filesep DTI_b0_edc_1_brain_mask ' -c 4 -p ' dir_DTI_pping filesep  'T2_temp_in_' DTI_sub_base '_b0_p%d.nii.gz -r 0.5 -w 0.05 -o ' ...
        dir_DTI_pping filesep DTI_sub_base '_edc1_b0_atropos_']); % here the b0 is the primary image but we can also do this with the T1 being the primary
    
    % % Dual channel Atropos segmentation of T1 & B0
    % % First bring your T1 template to diffusion space
    %         unix(['antsRegistrationSyNQuick.sh -d 3 -m ' T1_template_brain ...
    %         ' -f ' dir_DTI_pping filesep DTI_T1_6 ' -o ' ...
    %         dir_DTI_pping filesep 't1_template_2_B0_SyN_' DTI_sub_base  ' -x ' dir_DTI_pping filesep DTI_T1_6_mask ' -n 7 -j 1 -t s']);
    %     unix(['WarpImageMultiTransform 3 ' ttpm_csf ' ' dir_DTI_pping filesep ttpm_csf_in_subj ...
    %         ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    %     unix(['WarpImageMultiTransform 3 ' ttpm_wm ' ' dir_tpms_t1 filesep ttpm_wm_in_subj ...
    %         ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    %     unix(['WarpImageMultiTransform 3 ' ttpm_gm_cx ' ' dir_tpms_t1 filesep ttpm_gm_cx_in_subj ...
    %         ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    %     unix(['WarpImageMultiTransform 3 ' ttpm_gm_bg ' ' dir_tpms_t1 filesep ttpm_gm_bg_in_subj ...
    %         ' -R ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '1Warp.nii.gz ' ...
    %         dir_tpms_t1 filesep Temp_brain2subj_brain '0GenericAffine.mat --use-BSpline']);
    %     % T1_brain segmentation into Cortex, BG, WM,& CSF
    %     unix(['antsAtroposN4.sh -d 3 -a ' dir_brains_t1 filesep T1_subj_brain 'BrainExtractionBrain.nii.gz -x ' ...
    %         dir_brains_t1 filesep T1_subj_brain 'BrainExtractionMask.nii.gz -c 4 -p ' ...
    %         dir_tpms_t1 filesep (T1_base_name(1:end-12)) '%d.nii.gz -r 0.5 -w 0.05 -o ' dir_tpms_t1 filesep T1_base_name 'atropos_' ]);
    
    
    % Make an empty image for Vol5 to satisfy 5ttcheck
    unix(['fslmaths ' dir_DTI_pping filesep DTI_sub_base '_edc1_b0_atropos_SegmentationPosteriors3.nii.gz -mul 0 ' dir_DTI_pping filesep b0_empty_image ]);
    % Merge the 5 tpms
    unix(['fslmerge -t ' dir_DTI_pping filesep b0_pseudo_5ttgen ' ' dir_DTI_pping filesep b0_gm_cx ' ' dir_DTI_pping filesep b0_gm_bg ' ' ...
        dir_DTI_pping filesep b0_wm ' ' dir_DTI_pping filesep b0_csf ' ' dir_DTI_pping filesep b0_empty_image]);
    % 5tt check and 5ttgmwmi - to generate a gmwmi
    unix(['source ~/.bash_profile ; 5ttcheck -info -force -nthreads 7 -masks ' dir_DTI_pping filesep failed_5tt ' ' dir_DTI_pping filesep b0_pseudo_5ttgen]);
    unix(['source ~/.bash_profile ; 5tt2gmwmi -force -nthreads 7 ' dir_DTI_pping filesep  b0_pseudo_5ttgen ' ' dir_DTI_pping filesep b0_p5tt_gmwmi]);
    
    % will need to add smoothing here I think to do the GMWMI analysis part
    
    % All the stuff below is only commented out for testing the above steps
    % R 15/10/2018
    
    % This is aimed to go to fod and csd using the two shell multi
    % tissue CSD framework
    % before being able to do dwi2response I need a diffusion tissue
    % segmentation in distorted dwi space .... hmmm, I can do this with a
    % post=segmentation warp from the edc diffusion space, should be simple
    % enough... especially if I just use the 4D_5tt style segmentation.
    
    native_b0_p5tt_4D = ([DTI_sub_base '_dist_b0_p5tt_4D.nii.gz']);
    unix(['WarpTimeSeriesImageMultiTransform 4 ' dir_DTI_pping filesep b0_pseudo_5ttgen ' ' dir_DTI_pping filesep native_b0_p5tt_4D ...
        ' -R ' dir_DTI_pping filesep DTI_FA_5 ' -i ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN0GenericAffine.mat ' ...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN1InverseWarp.nii.gz --use-BSpline']);
    % do dhollander dwi2response
    dh_rf_wm = [DTI_sub_base '_Dh_RF_sfwm.txt'];
    dh_rf_gm = [DTI_sub_base '_Dh_RF_gm.txt'];
    dh_rf_csf = [DTI_sub_base '_Dh_RF_csf.txt'];
    dh_rf_voxels = [DTI_sub_base '_Dh_RF_voxels.mif'];
    csd_msmt_wm_fod = [DTI_sub_base '_msmt_wm_fod.mif'];
    csd_msmt_csf_fod = [DTI_sub_base '_msmt_csf_fod.mif'];
    wm_fod_id_warp = [DTI_sub_base 'wm_fod_id_warp.nii.gz'];
%     csf_fod_id_warp = [DTI_sub_base 'csf_fod_id_warp.nii.gz'];
    wm_fod_mrtrix_warp = [DTI_sub_base 'wm_fod_mrtrix_warp.nii.gz'];
    wm_fod_mrtrix_warp_corr = [DTI_sub_base 'wm_fod_mrtrix_warp_corr.nii.gz'];
    wm_fod_corrected = [DTI_sub_base 'wm_fod_edc_corrected.mif'];
    wm_fod_corr_mtn = [DTI_sub_base 'wm_fod_edc_corr_mtn.mif'];
%     csf_fod_mrtrix_warp = [DTI_sub_base 'csf_fod_mrtrix_warp.nii.gz'];
%     csf_fod_mrtrix_warp_corr = [DTI_sub_base 'csf_fod_mrtrix_warp_corr.nii.gz'];
%     csf_fod_corrected = [DTI_sub_base 'csf_fod_edc_corrected.mif'];
%     csf_fod_corr_mtn = [DTI_sub_base 'csf_fod_edc_corr_mtn.mif'];
    subj_wbrain_msmt_csd_tck_sfit1 = [DTI_sub_base '_wbrain_msmt_csd_tck_sfit1.tck'];
    
    unix(['source ~/.bash_profile ; dwi2response dhollander -force -nthreads 7 -voxels ' ...
        dir_DTI_pping filesep  dh_rf_voxels ' -mask ' dir_DTI_pping filesep DTI_BET_mask ' ' dir_DTI_pping filesep DTI_4 ' ' ...
        dir_DTI_pping filesep dh_rf_wm ' ' dir_DTI_pping filesep dh_rf_gm ' ' dir_DTI_pping filesep dh_rf_csf]);
    % do dwi2fod using wm and csf RF only
    unix(['source ~/.bash_profile ; dwi2fod -force -nthreads 7 msmt_csd -mask ' dir_DTI_pping filesep DTI_BET_mask ' '  dir_DTI_pping filesep DTI_4 ' ' ...
        dir_DTI_pping filesep dh_rf_wm  ' ' dir_DTI_pping filesep csd_msmt_wm_fod ' ' dir_DTI_pping filesep dh_rf_csf ' ' dir_DTI_pping filesep csd_msmt_csf_fod ]);
    % generate an mrtrix compatible warp from the ants warp and affine.mat
    % files to correct the fod and resulting tractogram for the effect of
    % SyN based EDC
    unix(['source ~/.bash_profile ; warpinit -force -nthreads 7 ' dir_DTI_pping filesep csd_msmt_wm_fod ' ' dir_DTI_pping filesep wm_fod_id_warp]);
    unix(['WarpTimeSeriesImageMultiTransform 4 ' dir_DTI_pping filesep wm_fod_id_warp ' ' dir_DTI_pping filesep wm_fod_mrtrix_warp ...
        ' -R ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyNWarped.nii.gz ' dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN1Warp.nii.gz ' ...
        dir_DTI_pping filesep DTI_sub_base '_FA_2_T1_SyN0GenericAffine.mat']);
    % Use these warps to rotate and correct the FODs for effect of SyN EDC
    unix(['source ~/.bash_profile ; warpcorrect -force -nthreads 7 '  dir_DTI_pping filesep wm_fod_mrtrix_warp ' ' dir_DTI_pping filesep wm_fod_mrtrix_warp_corr]);
    unix(['source ~/.bash_profile ; mrtransform -force -nthreads 7  -modulate -warp ' dir_DTI_pping filesep wm_fod_mrtrix_warp_corr ' -template ' dir_DTI_pping filesep DTI_T1_6 ' ' ...
        dir_DTI_pping filesep csd_msmt_wm_fod ' ' dir_DTI_pping filesep wm_fod_corrected ]); % my input warp needs to be a 4D one... hmf
    % Now we're ready for msmt_csd tractography with iFOD2
    unix(['source ~/.bash_profile ; mtnormalise -force -nthreads 7 -mask ' dir_DTI_pping filesep DTI_T1_6_mask ' ' ...
        dir_DTI_pping filesep wm_fod_corrected ' ' dir_DTI_pping filesep wm_fod_corr_mtn]);
    unix(['source ~/.bash_profile ; tckgen -force -algorithm iFOD2 -act ' dir_DTI_pping filesep b0_pseudo_5ttgen ' -seed_gmwmi ' ...
        dir_DTI_pping filesep b0_p5tt_gmwmi ' -mask ' dir_DTI_pping filesep DTI_T1_6_mask ' -select 2000000 -nthreads 7 '...
        dir_DTI_pping filesep wm_fod_corr_mtn ' ' dir_DTI_ready filesep subj_wbrain_msmt_csd_tck]); % works fine, but need to up streamlines no. I think...
    % SIFT the resulting tractogram
    unix(['source ~/.bash_profile ; tcksift -force -nthreads 7 -term_number 1000000 -act ' dir_DTI_pping filesep b0_pseudo_5ttgen ' '...
        dir_DTI_ready filesep subj_wbrain_msmt_csd_tck ' ' dir_DTI_pping filesep wm_fod_corr_mtn ' ' dir_DTI_ready filesep subj_wbrain_msmt_csd_tck_sfit1]);
    unix(['source ~/.bash_profile ; tcksift2 -force -nthreads 7 -act ' dir_DTI_pping filesep b0_pseudo_5ttgen ' ' dir_DTI_ready filesep subj_wbrain_msmt_csd_tck ...
        ' ' dir_DTI_pping filesep wm_fod_corr_mtn ' ' dir_DTI_ready filesep 'tract_weights.txt']);
    % Do we need SIFT1 as well ?
    % mtnormalise to correct for any remaining signal abnormalities or
    % biases
    % Use tck2connetome (weights, if possible FA, NOS, AFD, (FDI?), MD) -
    % maybe look at graph theory metrics also (rabbit hole! danger
    % danger!!)
    % need also fod2dec and fod2fixel -afd
    % connectome2tck find the fiber bundle of interest.
    
    
    %     % WON'T BE USING THIS FOR NOW.. MAYBE FOR A METHODS COMPARISON MANUSCRIPT ?... LATER
    %     % Probabilistic Tensor tractography using -act and -seedgmwmi with
    %     % 5 million streamlines
        unix(['source ~/.bash_profile ; tckgen -force -algorithm Tensor_Prob -act ' dir_DTI_pping filesep b0_pseudo_5ttgen ' -seed_gmwmi ' ...
        dir_DTI_pping filesep b0_p5tt_gmwmi ' -mask ' dir_DTI_pping filesep DTI_T1_6_mask ' -select 2000000 -nthreads 7 '...
        dir_DTI_pping filesep DTI_7_mif ' ' dir_DTI_ready filesep subj_wbrain_prob_tck]); % works fine, but need to up streamlines no. I think...
    %     % Bringing Atlas templates and labels to dMRI space
    %     unix(['antsRegistrationSyNQuick.sh -d 3 -m ' aal_template_brain ...
    %     ' -f ' dir_DTI_pping filesep DTI_T1_6 ' -o ' ...
    %     dir_DTI_pping filesep 'MNI1mbrain_2_' DTI_sub_base 'brain  -x ' ...
    %     dir_DTI_pping filesep DTI_T1_6_mask ' -n 7 -j 1 -t s']);
    %     unix(['antsRegistrationSyNQuick.sh -d 3 -m ' T1_template_brain ... % there was a typo here, now it should be okay 5/10/2018 R
    %     ' -f ' dir_DTI_pping filesep DTI_T1_6 ' -o ' ...
    %     dir_DTI_pping filesep 'PD25brain_2_' DTI_sub_base 'brain  -x ' ...
    %     dir_DTI_pping filesep DTI_T1_6_mask ' -n 7 -j 1 -t s']);
    %     unix(['WarpImageMultiTransform 3 ' aal_labels ' ' dir_DTI_ready filesep  aal_in_subj_DTI ...
    %     ' -R ' dir_DTI_pping filesep DTI_T1_6 ' ' ...
    %     dir_DTI_pping filesep 'MNI1mbrain_2_' DTI_sub_base 'brain1Warp.nii.gz ' dir_DTI_pping filesep 'MNI1mbrain_2_' DTI_sub_base 'brain0GenericAffine.mat --use-NN']); % still using NN as ML is rubbish
    %     unix(['WarpImageMultiTransform 3 ' BN_246_3D_labels ' ' dir_DTI_ready filesep BN246_in_subj_DTI ...
    %     ' -R ' dir_DTI_pping filesep DTI_T1_6 ' ' ...
    %     dir_DTI_pping filesep 'PD25brain_2_' DTI_sub_base 'brain1Warp.nii.gz ' dir_DTI_pping filesep 'PD25brain_2_' DTI_sub_base 'brain0GenericAffine.mat --use-NN']);
    %
    %     % fslsplit will be necessary before reorienting my tensors
    %     % This is what I need to import the DTs into ANTs
    %     %     ImageMath 4 dtiComp.nii.gz TimeSeriesDisassemble
    %     %     PT_001_DTI8_DT_edc1.nii.gz
    %     % i=0
    %     % for index in xx xy xz yy yz zz; do
    %     %    mv dtiComp100${i}.nii.gz dtiComp_${index}.nii.gz
    %     %    i=$((i+1))
    %     % done
    %     % ImageMath 3 dtAnts.nii.gz ComponentTo3DTensor dtiComp_
    %     % check https://sourceforge.net/p/advants/discussion/840261/thread/5e75969b/?limit=25
    %
    %     % b) then whole brain tractography
    %     % consider also try fsl and fslprobtrckx (maybe our data isn't good
    %     % enough)
    %     % c) use freesurfer segmentation for tck2connectome and connectome2tck
    %     % d) seeds of interest I imagine will include the hippo.s, the PhipG,
    %     % the PCC and PC, perhaps the supracalcarine cortex also ?
end




diary off

