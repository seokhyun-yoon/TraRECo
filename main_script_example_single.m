
%% Paired End Reads Assembly Example
% This script shows an assembly example of paired end read
% Input files
%   (1) Assembly Configuration: "config_mm9pe.par"
%   (2) Single-ended read: "Data/mm9_hydrolysis_c15pe_1f.fastq"

clear all;

%% Pre-compile C-functions (required only once)
compile_cmex;

%% Some input Arguments to TraRECo main
config_file = 'config_example_single.par';  % Assembly Config. file
Min_Tr_CvgDepth = [0 1 2 4 6 8];    % Read Coverage depth threshold used for Transcript selection 

%% TraRECo main
Out_file_prefix = trareco_v064( config_file, Min_Tr_CvgDepth );

%% Run BLAST-N 
% This requires BLAST Software installed on your computer
path_to_ref_transcriptome_db = 'sample_data/mm_c15_paired_ref.fasta';
Target_Cov = [0.8 0.9 0.95];   % Target coverages to measure
for k = 1:length(Min_Tr_CvgDepth)
    file_name = sprintf('%s_MinCvgDepth_%3.1f.fasta', Out_file_prefix, Min_Tr_CvgDepth(k) );
    run_blast( file_name, path_to_ref_transcriptome_db, Target_Cov );
end

