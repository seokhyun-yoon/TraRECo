
function cfg = trareco_get_config( fname_par, rd_length )

cfg.num_phs_loop = -1;
cfg.down_sample_factor = -1;
cfg.contig_ats = -1;
cfg.contig_ttl = -1;
cfg.max_pipeline = -1;         
cfg.num_pipeline = -1;         
cfg.max_contig_length = -1;
cfg.min_cd_ungrown = -1;
cfg.num_nmer_div = -1;
cfg.b_nmer_map_s = -1;
cfg.max_reads_to_load_m = -1;
cfg.bool_read_filter = -1;
cfg.bool_ss_ind = -1;

cfg.norm_dist_threshold = -1;
if exist('rd_length', 'var')
    cfg.nominal_read_length = rd_length;
else
    cfg.nominal_read_length = -1; 
    rd_length = -1;
end
cfg.min_overlap_depth = -1;
cfg.min_read_length = -1;
cfg.max_overlap_depth = -1;
cfg.min_tail_length = -1;
cfg.cntg_sel_threshold = -1;
cfg.num_nmer_div_js = -1;

cfg.junction_backoff = -1;
cfg.min_cvg_depth_js = -1;
cfg.min_cntg_length_js = -1;
cfg.safe_overlap_threshold = -1;
cfg.connection_threshold = -1;
cfg.connection_threshold_short = -1;
cfg.short_seg_threshold = -1;
cfg.max_n_cntg_per_group = -1;

cfg.min_cvg_depth = -1;
cfg.min_seg_length = -1;
cfg.max_num_paths = -1;
cfg.max_num_isoforms = -1;
cfg.b_tail_suppress = -1;
cfg.b_split_merged = -1;
cfg.mse_margin_percent = -1;
cfg.lasso_n_loops = -1;
cfg.lasso_step_size = -1;
cfg.min_tr_length = -1;
cfg.proc_unit_size = -1;
cfg.max_num_csets = -1;

cfg.read_mode = -1;
cfg.input_file_1 = -1;
cfg.input_file_2 = -1;
cfg.output_dir = '';
cfg.output_prefix = '';
cfg.memory_size_gb = -1;
cfg.min_rlen_percent = -1;
cfg.mean_frag_length = -1;

fname = fname_par;
fp = fopen( fname, 'r' );
if fp == -1
    fprintf('\nERROR: Cannot open User Cofigutaion. Check your configuration file %s ', fname );
    cfg = -1;
else
    while(1)
       str = fgets(fp);
       if str < 0
           break;
       else
           if length(str) > 5 
               kword = 'SYSTEM_MEMORY_SIZE_GB';
               if strcmp(str(1:length(kword)), kword)
                   cfg.memory_size_gb = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'OUTPUT_DIRECTORY';
               if strcmp(str(1:length(kword)), kword)
                   cfg.output_dir = sscanf(str(length(kword)+1:end), '%s');
               end
               kword = 'OUTPUT_FILE_PREFIX';
               if strcmp(str(1:length(kword)), kword)
                   cfg.output_prefix = sscanf(str(length(kword)+1:end), '%s');
               end
               kword = 'INPUT_FILE_SINGLE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.input_file_1 = sscanf(str(length(kword)+1:end), '%s');
                   cfg.read_mode = cfg.read_mode + 1;
               end
               kword = 'INPUT_FILE_PAIRED_1';
               if strcmp(str(1:length(kword)), kword)
                   cfg.input_file_1 = sscanf(str(length(kword)+1:end), '%s');
                   cfg.read_mode = cfg.read_mode + 1;
               end
               kword = 'INPUT_FILE_PAIRED_2';
               if strcmp(str(1:length(kword)), kword)
                   cfg.input_file_2 = sscanf(str(length(kword)+1:end), '%s');
                   cfg.read_mode = cfg.read_mode + 1;
               end
               kword = 'NUM_PHS_LOOPS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.num_phs_loop = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_NUM_READS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.num_reads_total = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'DOWN_SAMPLE_FACTOR';
               if strcmp(str(1:length(kword)), kword)
                   cfg.down_sample_factor = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'CONTIG_TIME_TO_LIVE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.contig_ttl = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'CG_WINDOW_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.contig_ats = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_CD_TO_REMOVE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cd_ungrown = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MAX_NUM_PIPELINES';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_pipeline = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_READ_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_read_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_TAIL_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_tail_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'NOMINAL_READ_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.nominal_read_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_CONTIG_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_contig_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_OVERLAP_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_overlap_depth = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'NORM_DIST_THRESHOLD';
               if strcmp(str(1:length(kword)), kword)
                   cfg.norm_dist_threshold = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'NMER_DOWN_SAMPLE_CG';
               if strcmp(str(1:length(kword)), kword)
                   cfg.num_nmer_div = sscanf(str(length(kword)+1:end), '%d');
               end
               
               kword = 'STRAND_SPECIFICITY';
               if strcmp(str(1:length(kword)), kword)
                   cfg.bool_ss_ind = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'READ_FILTER_ENABLE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.bool_read_filter = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_RDS_TO_LOAD_M';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_reads_to_load_m = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'NMER_DOWN_SAMPLE_JS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.num_nmer_div_js = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'B_NMER_MAP_TYPE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.b_nmer_map_s = sscanf(str(length(kword)+1:end), '%d');
               end
                kword = 'MIN_SEL_THRESHOLD';
               if strcmp(str(1:length(kword)), kword)
                   cfg.cntg_sel_threshold = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'JUNCTION_BACKOFF';
               if strcmp(str(1:length(kword)), kword)
                   cfg.junction_backoff = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'JUNCTION_THRESHOLD';
               if strcmp(str(1:length(kword)), kword)
                   cfg.connection_threshold = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'CONTIG_LENGTH_SHORT';
               if strcmp(str(1:length(kword)), kword)
                   cfg.short_seg_threshold = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'S_OVERLAP_THRESHOLD';
               if strcmp(str(1:length(kword)), kword)
                   cfg.safe_overlap_threshold = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'CONNECTION_THRESHOLD';
               if strcmp(str(1:length(kword)), kword)
                   cfg.connection_threshold_short = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_COVERAGE_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cvg_depth_js = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_CONTIG_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cvg_depth_js = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_CONTIG_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cntg_length_js = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_TRANSCRIPT_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cvg_depth = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_TR_CAND_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cvg_depth = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_SEGMENT_DEPTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_cvg_depth = sscanf(str(length(kword)+1:end), '%f');
               end
               kword = 'MIN_SEGMENT_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_seg_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_N_CNTGS_PER_GROUP';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_n_cntg_per_group = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MSE_MARGIN_PERCENT';
               if strcmp(str(1:length(kword)), kword)
                   cfg.mse_margin_percent = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'LASSO_NUM_OF_LOOPS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.lasso_n_loops = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'BOOL_TAIL_SUPPRESS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.b_tail_suppress = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'BOOL_SPLIT_MERGED';
               if strcmp(str(1:length(kword)), kword)
                   cfg.b_split_merged = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_NUM_PATHS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_num_paths = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_NUM_CAND_SETS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_num_csets = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MAX_NUM_ISOFORMS';
               if strcmp(str(1:length(kword)), kword)
                   cfg.max_num_isoforms = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_TR_CAND_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_tr_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_TRANSCRIPT_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_tr_length = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'PROCESS_UNIT_SIZE';
               if strcmp(str(1:length(kword)), kword)
                   cfg.proc_unit_size = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MIN_TLEN_PERCENT_WRT_MAX';
               if strcmp(str(1:length(kword)), kword)
                   cfg.min_rlen_percent = sscanf(str(length(kword)+1:end), '%d');
               end
               kword = 'MEAN_FRAGMENT_LENGTH';
               if strcmp(str(1:length(kword)), kword)
                   cfg.mean_frag_length = sscanf(str(length(kword)+1:end), '%d');
               end
           end
       end
    end
    fclose(fp);
    
    if cfg.nominal_read_length <= 0 
        fprintf('\nERROR: NOMINAL_READ_LENGTH not specified in your configuration file ' );
        fprintf('\nPlease specify nominal read length with the keyworrd ''NOMINAL_READ_LENGTH''' );
        cfg = -1;
        return;
    else
        if isempty(cfg.output_dir)
            cfg.output_dir = 'trareco';
        end
        if isempty(cfg.output_prefix)
            cfg.output_prefix = 'trareco';
        end
        if cfg.nominal_read_length < 0 && rd_length > 0
            cfg.nominal_read_length = rd_length;
        else
            if cfg.nominal_read_length > 0 && rd_length < 0
                % rd_length = cfg.nominal_read_length;
            else
                if cfg.nominal_read_length ~= rd_length
                    fprintf('\nWARNING: System detected nominal read length = %d ', rd_length );
                    fprintf('\nWhile, you specified (in your configuration) NOMINAL_READ_LENGTH = %d ', cfg.nominal_read_length );
                end
            end
        end
    end
    if cfg.norm_dist_threshold < 0
        cfg = get_default_values( cfg.nominal_read_length );
        return;
    else
        cfg_default = get_default_values( cfg.nominal_read_length, cfg.norm_dist_threshold );
    end
    if cfg.num_phs_loop < 0
        cfg.num_phs_loop = cfg_default.num_phs_loop;
    end
    if cfg.down_sample_factor < 0
        cfg.down_sample_factor = cfg_default.down_sample_factor;
    end
    if cfg.min_cd_ungrown < 0
        cfg.min_cd_ungrown = cfg_default.min_cd_ungrown;
    end
    if cfg.bool_read_filter < 0
        cfg.bool_read_filter = cfg_default.bool_read_filter;
    end
    if cfg.bool_ss_ind < 0
        cfg.bool_ss_ind = cfg_default.bool_ss_ind;
    end
    if cfg.contig_ats < 0
        cfg.contig_ats = cfg_default.contig_ats;
    end
    if cfg.contig_ttl < 0
        cfg.contig_ttl = cfg_default.contig_ttl;
    end
    if cfg.max_pipeline < 0
        cfg.max_pipeline = cfg_default.max_pipeline;  
    end
    if cfg.num_pipeline < 0
        cfg.num_pipeline = cfg_default.num_pipeline;   
    end
    if cfg.proc_unit_size < 0
        cfg.proc_unit_size = cfg_default.proc_unit_size;
    end
    if cfg.memory_size_gb < 0
        cfg.memory_size_gb = cfg_default.memory_size_gb; 
    else
        if cfg.memory_size_gb < 16
            cfg.max_pipeline = round( cfg.max_pipeline*cfg.memory_size_gb/16 );
            cfg.num_pipeline = round( cfg.num_pipeline*cfg.memory_size_gb/16 );
            cfg.proc_unit_size = round( cfg.proc_unit_size*cfg.memory_size_gb/16 );
        end
    end
    if cfg.max_contig_length < 0
        cfg.max_contig_length = cfg_default.max_contig_length;
    end
    if cfg.min_overlap_depth < 0
        cfg.min_overlap_depth = cfg_default.min_overlap_depth;
    end
    if cfg.max_overlap_depth < 0
        cfg.max_overlap_depth = cfg_default.max_overlap_depth;
    end
    if cfg.min_read_length < 0
        cfg.min_read_length = cfg_default.min_read_length;
    end
    if cfg.min_tail_length < 0
        cfg.min_tail_length = cfg_default.min_tail_length;
    end
    if cfg.num_nmer_div < 0
        cfg.num_nmer_div = cfg_default.num_nmer_div;
    end
    if cfg.max_reads_to_load_m < 0
        cfg.max_reads_to_load_m = cfg_default.max_reads_to_load_m;
    end
    if cfg.num_nmer_div_js < 0
        cfg.num_nmer_div_js = cfg_default.num_nmer_div_js;
    end
    if cfg.b_nmer_map_s < 0
        cfg.b_nmer_map_s = cfg_default.b_nmer_map_s;
    end
    if cfg.junction_backoff < 0
        cfg.junction_backoff = cfg_default.junction_backoff;
    end
    if cfg.min_cvg_depth_js < 0
        cfg.min_cvg_depth_js = cfg_default.min_cvg_depth_js;
    end
    if cfg.min_cntg_length_js < 0
        cfg.min_cntg_length_js = cfg_default.min_cntg_length_js;
    end
    if cfg.safe_overlap_threshold < 0
        cfg.safe_overlap_threshold = cfg_default.safe_overlap_threshold;
    end
    if cfg.connection_threshold < 0
        cfg.connection_threshold = cfg_default.connection_threshold;
    end
    if cfg.connection_threshold_short < 0
        cfg.connection_threshold_short = cfg_default.connection_threshold_short;
    end
    if cfg.short_seg_threshold < 0
        cfg.short_seg_threshold = cfg_default.short_seg_threshold;
    end
    if cfg.min_tail_length < 0
        cfg.min_tail_length = cfg_default.min_tail_length;
    end
    if cfg.cntg_sel_threshold < 0
        cfg.cntg_sel_threshold = cfg_default.cntg_sel_threshold;
    end
    if cfg.min_cvg_depth < 0
        cfg.min_cvg_depth = cfg_default.min_cvg_depth;
    end
    if cfg.min_seg_length < 0
        cfg.min_seg_length = cfg_default.min_seg_length;
    end
    if cfg.max_n_cntg_per_group < 0
        cfg.max_n_cntg_per_group = cfg_default.max_n_cntg_per_group;
    end
    if cfg.b_tail_suppress < 0
        cfg.b_tail_suppress = cfg_default.b_tail_suppress;
    end
    if cfg.b_split_merged < 0
        cfg.b_split_merged = cfg_default.b_split_merged;
    end
    if cfg.mse_margin_percent < 0
        cfg.mse_margin_percent = cfg_default.mse_margin_percent;
    end
    if cfg.lasso_n_loops < 0
        cfg.lasso_n_loops = cfg_default.lasso_n_loops;
    end
    if cfg.lasso_step_size < 0
        cfg.lasso_step_size = cfg_default.lasso_step_size;
    end
    if cfg.max_num_paths < 0
        cfg.max_num_paths = cfg_default.max_num_paths;
    end
    if cfg.max_num_isoforms < 0
        cfg.max_num_isoforms = cfg_default.max_num_isoforms;
    end
    if cfg.min_tr_length < 0
        cfg.min_tr_length = cfg_default.min_tr_length;
    end
    if cfg.max_num_csets < 0
        cfg.max_num_csets = cfg_default.max_num_csets;
    end
    if cfg.min_rlen_percent < 0
        cfg.min_rlen_percent = cfg_default.min_rlen_percent;
    end
    if cfg.mean_frag_length < 0
        if cfg.read_mode == 0
            cfg.mean_frag_length = cfg.nominal_read_length;
        else
            cfg.mean_frag_length = cfg.nominal_read_length*3;
        end
    end
end
end

function cfg = get_default_values( rd_len, dist_th )

    cfg.num_phs_loop = 120;
    cfg.down_sample_factor = 1;
    cfg.contig_ats = 0;
    cfg.contig_ttl = 300000;
    cfg.max_pipeline = 60000;  
    cfg.num_pipeline = 60000;   
    cfg.proc_unit_size = 15000;
    cfg.max_contig_length = 3000;
    if exist('dist_th', 'var')
        cfg.norm_dist_threshold = dist_th;
    else
        cfg.norm_dist_threshold = 0.06;
    end
    cfg.min_cd_ungrown = 2.1;
    cfg.max_reads_to_load_m = 0;
    cfg.nominal_read_length = rd_len;
    cfg.min_overlap_depth = 44; 
    cfg.max_overlap_depth = cfg.nominal_read_length;
    cfg.min_read_length = round(cfg.nominal_read_length*0.6);
    cfg.num_nmer_div = 8;
    cfg.b_nmer_map_s = 0;
    cfg.bool_read_filter = 1;
    cfg.bool_ss_ind = 0;
    cfg.min_cvg_depth_js = 1.8;
    cfg.min_cntg_length_js = round(cfg.nominal_read_length*0.8);
    cfg.safe_overlap_threshold = 50;
    cfg.connection_threshold = 24; 
    cfg.connection_threshold_short = 24; 
    cfg.short_seg_threshold = round(cfg.nominal_read_length);
    cfg.max_n_cntg_per_group = 40;
    cfg.num_nmer_div_js = 4;
    cfg.min_tail_length = min( max(22, round(cfg.connection_threshold/2) ), cfg.connection_threshold_short );
    cfg.cntg_sel_threshold = 2;
    cfg.junction_backoff = 12;
    cfg.min_cvg_depth = 0;
    cfg.min_seg_length = cfg.nominal_read_length - cfg.connection_threshold;
    cfg.max_num_paths = 150;
    cfg.max_num_isoforms = 200;
    cfg.b_tail_suppress = 1;
    cfg.b_split_merged = 0;
    cfg.mse_margin_percent = 3;
    cfg.lasso_n_loops = 400;
    cfg.lasso_step_size = 0.5/cfg.lasso_n_loops;
    cfg.min_tr_length = 200;
    cfg.max_num_csets = 1;
    cfg.memory_size_gb = 16;
    cfg.min_rlen_percent = 40;
end

