
%% Paired End Reads Assembly Example
% This script shows an assembly example of paired end read
% Input files
%   (1) Assembly Configuration: "config_mm9pe.par"
%   (2) Paired read 1: "mm9_hydrolysis_c15pe_1f.fastq"
%   (3) Paired read 2: "mm9_hydrolysis_c15pe_2r.fastq"

clear all;

%% Pre-compile C-functions (required only once)
compile_cmex;

%% Some input Arguments to TraRECo main
config_file = 'config_example_paired.par';  % Assembly Config. file
Min_Tr_CvgDepth = [0 1 2 4 6 8 12 16];    % Read Coverage depth threshold used for Transcript selection 
Min_Tr_Length = 200;
mode = 1;

%% TraRECo main
Out_file_prefix = trareco_v067( config_file, Min_Tr_CvgDepth, mode );

%% Run BLAST-N 
% This requires BLAST Software installed on your computer
path_to_ref_transcriptome_db = 'sample_data/mm_c15_paired_ref.fasta';
Target_Cov = [0.8 0.9 0.95];   % Target coverages to measure
for k = 1:length(Min_Tr_CvgDepth)
    file_name = sprintf('%s_MinCvgDepth_%d.fasta', Out_file_prefix, Min_Tr_CvgDepth(k) );
    tr_cands = run_blast( file_name, path_to_ref_transcriptome_db, Target_Cov, Min_Tr_Length );
    
    abn_rpkm_true = [tr_cands(1:end).abn_rpkm_true]';
    abn_rpkm_est = [tr_cands(1:end).abn_rpkm_est]';
    [~, gof] = fit(abn_rpkm_true, abn_rpkm_est,'poly1');
    R2 = gof.rsquare;
    str = sprintf('R^2 = %f', R2 );
    disp(str);
end

