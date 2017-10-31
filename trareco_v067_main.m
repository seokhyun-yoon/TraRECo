
function trareco_v067_main( an1, arg1, an2, arg2, an3, arg3, an4, arg4 )

Min_Tr_Length = 0;
path_to_ref_transcriptome_db = [];
config_file = [];
blast_flag = 0;
run_mode = 1;
if ~exist( 'an1', 'var' ) || ~exist( 'arg1', 'var' )
    n_args = 0;
else
    n_args = 1;
    if ~exist( 'an2', 'var' ) || ~exist( 'arg2', 'var' )
    else
        n_args = n_args + 1;
        if ~exist( 'an3', 'var' ) || ~exist( 'arg3', 'var' )
        else
            n_args = n_args + 1;
            if ~exist( 'an4', 'var' ) || ~exist( 'arg4', 'var' )
            else
                n_args = n_args + 1;
            end
        end
    end
end

if n_args < 1
    fprintf('   Not enough arguments provided ... \n');
    print_usage( )
    return;
else
    for n = 1:n_args
        switch( n )
            case 1, 
                an = an1;
                arg = arg1;
            case 2,
                an = an2;
                arg = arg2;
            case 3,
                an = an3;
                arg = arg3;
            case 4,
                an = an4;
                arg = arg4;
            otherwise, break;
        end
        at = get_option( an );
        switch( at )
            case 1, 
                config_file = arg;
            case 2,
                path_to_ref_transcriptome_db = arg;
            case 3,
                Min_Tr_Length = str2num( arg );
            case 4,
                run_mode = str2num( arg );
                if run_mode > 3
                    run_mode = 3;
                end
            otherwise
        end
    end
end
Min_Tr_CvgDepth = [0 1 2 4 6 8 12 16];  % Read Coverage depth threshold used for Transcript selection 
Target_Cov = [0.8 0.9 0.95];   % Target coverages to measure

if ~exist( config_file, 'file' ) 
    if ~exist( config_file, 'file' )
        fprintf('   Configuration file %d does not exist .. \n', config_file );
    end
else
    fprintf('Running TraRECo v0.67 .. \nConfiguration file: %s \n', config_file );
    if exist( path_to_ref_transcriptome_db, 'file' )
        fprintf('Ref. Transcriptome to compare: %s \n', path_to_ref_transcriptome_db );
        blast_flag = 1;
        fprintf('Run BLAST option set \n' );
        fprintf('BLASTN must be installed on your computer. \n' );
    end
    if run_mode > 1
        fprintf('Jump to Step: %d \n', run_mode );
    end
  
    %% TraRECo main
    Out_file_prefix = trareco_v067a( config_file, Min_Tr_CvgDepth, run_mode );
    %% Run BLAST-N 
    % This requires BLAST Software installed on your computer
    
    if blast_flag > 0
        for k = 1:length(Min_Tr_CvgDepth)
            file_name = sprintf('%s_MinCvgDepth_%d.fasta', Out_file_prefix, Min_Tr_CvgDepth(k) );
            tr_cands = run_blast( file_name, path_to_ref_transcriptome_db, Target_Cov, Min_Tr_Length );

            abn_rpkm_true = [tr_cands(1:end).abn_rpkm_true]';
            abn_rpkm_est = [tr_cands(1:end).abn_rpkm_est]';
            [~, gof] = fit(abn_rpkm_true, abn_rpkm_est,'poly1');
            R2 = gof.rsquare;
            str = sprintf('Abundance Estimation R^2 measure = %f', R2 );
            disp(str);
        end
    end
end

end

function out_file_name = trareco_v067a( config_file, Min_tr_abundance, mode )

%% TraRECo main
fprintf('+-----------------------------------+\n');
fprintf('|        TraRECo version 0.67       |\n');
fprintf('+-----------------------------------+\n');
if mode < 2
    trareco_cg_v067( config_file );
end
if mode < 3
    trareco_js_v067( config_file );  
end
out_file_name = trareco_td_v067( config_file );
out_file_name = trareco_tr_filter( out_file_name, Min_tr_abundance );
fprintf('+-----------------------------------+\n');
fprintf('|    Thank you for using TraRECo    |\n');
fprintf('+-----------------------------------+\n');

end

function print_usage( )
    fprintf('   TraRECo options ... \n');
    fprintf('   -c [path/configuration file name] \n');
    fprintf('   -r [path/reference transcriptome file name (in fasta format)] -- optional \n');
    fprintf('      To use -r option, standalone version of BLASTN must be installed on your computer \n');
end

function at = get_option( an1 )
    at = 0;
    if isempty( an1 )
        fprintf('   Wrong usage ... \n');
        print_usage( )
        return;
    else
        if an1(1) ~= '-'
            fprintf('   Wrong usage ... \n');
            print_usage( )
            return;
        else
            switch( an1(2) )
                case 'c', at = 1;
                case 'r', at = 2;
                case 'm', at = 3;
                case 'j', at = 4;
                otherwise, at = 0;
            end
        end
    end
end

