function fname_out = trareco_cg_v067( fname_param )

fprintf('Checking input files ... ');
cfg = trareco_get_config( fname_param );
if isnumeric(cfg) || cfg.read_mode < 0
    return;
end
fname_ext = cfg.input_file_1;
[N_reads_max, read_type, fname, fext, read_length] = est_n_reads( fname_ext );
if read_type == 0 || N_reads_max < 0
    fprintf('File type not supported or cannot open \n');
    return;
end
fname_to_save = cfg.output_prefix; 
if cfg.read_mode == 0
    fnamep = -1;
else
    fnamep_ext = cfg.input_file_2;
    [N_reads_maxp, read_typep, fnamep, fextp, read_lengthp] = est_n_reads( fnamep_ext );
    if read_typep ~= read_type || ~strcmp( fext, fextp )
        fprintf('Types of the two input files are different \n');
        return;
    end
    if read_lengthp ~= read_length || N_reads_maxp ~= N_reads_max
        % No problem for now
    end
end
% fprintf('done');

if read_type == 2
    str_type = 'FASTQ';
    f_type = 'q';
else
    if read_type == 1
        str_type = 'FASTA';
        f_type = 'a';
    else
        fprintf('File type not supported. \nSupported types: .fastq, .fq, .fasta, .fa \n');
        return;
    end
end

if cfg.nominal_read_length ~= read_length
    fprintf('\n   WARNING: System detected nominal read length = %d ', read_length );
    fprintf('\n   While, you specified (in your configuration) NOMINAL_READ_LENGTH = %d ', cfg.nominal_read_length );
end

% fname_fax = sprintf('%s.%s', fname, fext );
if length(fnamep) == 1 && fnamep(1) < 0
    N_reads_max = ceil(N_reads_max/100000)*100000;
    fprintf( '\n   Type: %s, Single, Approx. # of reads: %dK, Read length: %d \n', ...
        str_type, ceil(N_reads_max/1000), cfg.nominal_read_length );
else
    N_reads_max = ceil(N_reads_max/100000)*100000*2;
    fprintf( '\n   Type: %s, Paired, Approx. # of reads: %dK (x2), Read length: %d \n', ...
        str_type, ceil(N_reads_max/2000), cfg.nominal_read_length );
end

%% Simulation Parameters
sim_prm.N_reads_max = N_reads_max; %0;
sim_prm.down_sample_factor = cfg.down_sample_factor;
sim_prm.num_N_mer = 8;
if cfg.num_nmer_div >= 1 && cfg.num_nmer_div <= sim_prm.num_N_mer
    sim_prm.seg_sel_dsr = round(cfg.num_nmer_div);
else
    sim_prm.seg_sel_dsr = sim_prm.num_N_mer;
end
sim_prm.disp_period = 10000;
sim_prm.update_period = 10101;
sim_prm.b_inter_disp = 0;
Seg_Nmer_map = repmat( zeros( cfg.max_pipeline,1, 'uint8' ), 1, 4^sim_prm.num_N_mer );
if cfg.b_nmer_map_s > 0
    Seg_Nmer_map2 = Seg_Nmer_map;
    slen = cfg.nominal_read_length - sim_prm.num_N_mer + 1;
end
nr_time_start = 200000;

N_phs_loop = cfg.num_phs_loop;
Min_CV_to_remove = cfg.min_cd_ungrown; 
Q_ref = 2; 
Normalized_Dist_threshold = cfg.norm_dist_threshold; 
Max_pipeline = cfg.max_pipeline; 
Min_tail_length = cfg.min_tail_length;
Max_seg_length_for_combining = round(cfg.max_contig_length*2.25/5);
Down_sample_factor = sim_prm.down_sample_factor;
Seg_sel_dsr = sim_prm.seg_sel_dsr;
Min_ovlp_dep = cfg.min_overlap_depth;
cfg_val = [cfg.norm_dist_threshold Min_ovlp_dep cfg.max_overlap_depth cfg.bool_ss_ind];
Initial_buf_length = uint32(2000); 

%% Running Variables
seg_seq_tmp1 = int8( zeros(1,cfg.max_contig_length*4) );
seg_cvg_tmp1 = uint32( zeros(4,cfg.max_contig_length*4) );
seg_seq_tmp2 = int8( zeros(1,cfg.max_contig_length*4) );
seg_cvg_tmp2 = uint32( zeros(4,cfg.max_contig_length*4) );
mxv = zeros(1,cfg.max_contig_length*4, 'uint32');
mxi = zeros(1,cfg.max_contig_length*4, 'uint32');

sg.id = uint32(0);
sg.blen = uint32(Initial_buf_length);
sg.len = uint32(0);
sg.lin = uint32(0);
sg.stm = int32(0);
sg.seq = zeros(1,Initial_buf_length, 'int8');
sg.cvg_dep = zeros(4,Initial_buf_length,'uint32');
sg.ave_cvg_dep = uint32(0);
seg = repmat( sg, cfg.max_pipeline, 1 );

b_seg_valid = zeros(1,cfg.max_pipeline, 'int16');
priority = zeros(cfg.max_pipeline,1,'int16');
n_em_po_pdf = zeros(4,4);

if ~isempty( cfg.output_dir ) && ischar( cfg.output_dir )
    if exist( cfg.output_dir, 'dir' )
        fprintf( '   WARNING: The directory you specified already exists\n' );
        fprintf( '   The existing files will be overwritten.\n' );
    else
        mkdir( cfg.output_dir );
    end
    fname_out = sprintf('%s/%s', cfg.output_dir, fname_to_save ); 
else
    fname_out = sprintf('%s', fname_to_save ); 
end
fname_log = sprintf('%s.log', fname_out );
fp_log = fopen( fname_log, 'wt' );
if fp_log < 0
    fprintf('\nCannot open: %s \n', fname_log );
    return;
end

if cfg.read_mode == 0
    fprintf(fp_log, '# Input: %s.%s\n', fname, fext );
else
    fprintf(fp_log, '# Input 1: %s.%s\n', fname, fext );
    fprintf(fp_log, '# Input 2: %s.%s\n', fnamep, fext );
end
fprintf(fp_log, '# Assembly parameters: %s\n', fname_param );
fprintf(fp_log, 'Number of parallel contigs: %d\n', cfg.max_pipeline );
fprintf(fp_log, 'Nominal read length: %d\n', cfg.nominal_read_length );
fprintf(fp_log, 'Minimum overlap width: %d\n', Min_ovlp_dep );
fprintf(fp_log, 'Minimum read length: %d\n', cfg.min_read_length );
fprintf(fp_log, 'Normalized Distance threshold: %d\n', cfg.norm_dist_threshold );
if cfg.b_tail_suppress > 0
    fprintf(fp_log, 'Multiple Poly(A) sites suppressed \n' );
else
end

%% Build segment library

min_sel_th = cfg.cntg_sel_threshold;
Ave_num_error_normalized = ceil(Min_ovlp_dep*Normalized_Dist_threshold);
Seg_select_threshold = round( max( (Min_ovlp_dep - (sim_prm.num_N_mer-Seg_sel_dsr))/Seg_sel_dsr - Ave_num_error_normalized*(sim_prm.num_N_mer/Seg_sel_dsr), min_sel_th ) );
Ave_num_error_normalized_2 = ceil(cfg.nominal_read_length*Normalized_Dist_threshold);
Seg_select_threshold_2 = round( max( (cfg.nominal_read_length - (sim_prm.num_N_mer-Seg_sel_dsr))/Seg_sel_dsr - Ave_num_error_normalized_2*(sim_prm.num_N_mer/Seg_sel_dsr), min_sel_th ) );

dstr = char( datestr(now) );
if length(fnamep) == 1 && fnamep(1) < 0
    str_disp = sprintf('Contig growing started with CST(%d,%d,%d) ...  %s ', ...
        Seg_select_threshold, Seg_select_threshold_2, Min_ovlp_dep, dstr );
    R_mode = 0;
else
    str_disp = sprintf('Contig growing started with CST(%d,%d,%d) ...  %s ', ...
        Seg_select_threshold, Seg_select_threshold_2, Min_ovlp_dep, dstr );
    R_mode = 1;
end
disp(str_disp);
fprintf( fp_log, '\n%s\n', str_disp ); 

n_seg_found = 0;
n_seg_maxid = 0;
k_disp_flag = 0;
N_reads_valid = 0;
N_reads_processed = 0;
N_discarded = 0;

fname_txt = sprintf('%s.cntg', fname_out );
fp_t = fopen( fname_txt, 'wt' );
if fp_t < 0
    fprintf('\nFile not exists or cannot open: %s \n', fname_txt );
    return;
end
% State variables
N_seg = 0;
k = 0; 
k_first_drop = 0;
k_drop_time = 0;
Mode = 0;
Mode_Rd = 0;

n_base = 0;
sum_sum = 0;
sum_max = 0;

cntg_cnt = 0;
cntg_valid_ind = (-1).*ones( sim_prm.N_reads_max, 1 ); % -1: uninitialized
read_cntg_map = zeros( sim_prm.N_reads_max, 3 );

% sim_prm.N_reads_max = 100000;
for phs = 1:1:N_phs_loop

C_TTL = cfg.contig_ttl*phs;
if phs == 1
    fname_fa = sprintf('%s', fname_ext );
    fp_fa = fopen( fname_fa, 'r' );
    if fp_fa < 0
        fprintf('\nFile not exist or cannot open: %s \n', fname_fa );
        return;
    end
    if R_mode == 0
%         fname_s1 = sprintf('%s.reads', fname_out );
%         fp_s(1) = fopen( fname_s1, 'wt' );
    else
%         fname_s1 = sprintf('%s_1.reads', fname_out );
%         fp_s(1) = fopen( fname_s1, 'wt' );
%         fname_s2 = sprintf('%s_2.reads', fname_out );
%         fp_s(2) = fopen( fname_s2, 'wt' );
        
        fname_fp = sprintf('%s', fnamep_ext );
        fp_fp = fopen( fname_fp, 'r' );
        if fp_fp < 0
            fprintf('\nFile not exist or cannot open: %s \n', fname_fp );
            return;
        end
    end
    fname_fd = sprintf('%s_tmp%d.fasta', fname_out, phs );
    fp_fd = fopen( fname_fd, 'w' );
    if fp_fd < 0
        fprintf('\nCannot open: %s \n', fname_fd );
        return;
    end
    N_reads_max = sim_prm.N_reads_max;
    N_dropped = 0;
    DSF = Down_sample_factor;
else
    if phs < 2
        Min_ovlp_dep = cfg.min_overlap_depth;
    else
        if N_dropped*0.5 > C_TTL && Drop_rate < 0.9
            Min_ovlp_dep = min( cfg.min_overlap_depth, 44 );
        else
            max_cntg_len = max( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] );
            if max_cntg_len > cfg.nominal_read_length*5 && N_dropped*0.8 > C_TTL % cfg.nominal_read_length*5
                Min_ovlp_dep = min( cfg.min_overlap_depth, 36 );
            else
                Min_ovlp_dep = min( cfg.min_overlap_depth, 28 );
            end
        end
    end
    
    if Mode_Rd == 0
        if N_dropped >= 0 % cfg.max_reads_to_load_m*1e6
            fname_fa = sprintf('%s_tmp%d.fasta', fname_out, phs-1 );
            fp_fa = fopen( fname_fa, 'r' );
            if fp_fa < 0
                fprintf('\nFile not exist or cannot open: %s \n', fname_fa );
                return;
            end
            fname_fd = sprintf('%s_tmp%d.fasta', fname_out, phs );
            fp_fd = fopen( fname_fd, 'w' );
            if fp_fd < 0
                fprintf('\nCannot open: %s \n', fname_fd );
                return;
            end
        else
            Reads_tmp.r_or_c = 0;
            Reads_tmp.prev_id = 0;
            Reads_tmp.pei = 0;
            Reads_tmp.len = 0;
            Reads_tmp.seq = int8(0); %zeros(1,round(cfg.nominal_read_length*1.25),'int8');
            Reads_tmp.dep = 0;
            Reads_tmp.k_drop = 0;
            Reads_all = repmat( Reads_tmp, N_dropped, 1 );
            Dropped = repmat( Reads_tmp, N_dropped, 1 );
            
            fname_fa = sprintf('%s_tmp%d.fasta', fname_out, phs-1 );
            fp_fa = fopen( fname_fa, 'rt' );
            if fp_fa < 0
                fprintf('\nFile not exist or cannot open: %s \n', fname_fa );
                return;
            end
            fprintf('   Loading %s ... ', fname_fa );
            n_step = 1000; %round(N_dropped/20);
            Nchar = 0;
            for m = 1:1:N_dropped
                [rd_seq, rd_len, Cdep, r_or_c, prev_id, pei, k_drop_time] = f04_fasta_get_read_v09( fp_fa, cfg.min_read_length, DSF, 1, 0, Normalized_Dist_threshold, cfg.nominal_read_length, 0 );
                if rd_len > 0
                    Reads_all(m).r_or_c = r_or_c;
                    Reads_all(m).prev_id = prev_id;
                    Reads_all(m).pei = pei;
                    Reads_all(m).len = rd_len;
                    Reads_all(m).seq = rd_seq;
                    Reads_all(m).dep = Cdep;
                    Reads_all(m).k_drop = k_drop_time;
                    if mod( m, n_step ) == 0
                        if Nchar > 0
                            fprintf(repmat('\b', 1, Nchar));
                        end
                        Nchar = fprintf('%d/%d', m, N_dropped );
                    end
                else
                    N_dropped = m-1;
                    break;
                end
            end
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            fprintf('%d/%d', m, N_dropped );
            fprintf(' done \n' );
            fclose(fp_fa);
            Mode_Rd = 1;
        end
    else
        Reads_all = Dropped(1:N_dropped);
        Dropped = Reads_all;
    end
    
    N_reads_max = N_dropped;
    N_dropped = 0;
    DSF = 1;
end

kclk1 = clock();
xcntr = 0;
N_reads_sofar = 0;
pr_n = 0;
pr_g = 0;
n_seg_selected = 0;
n_seg_selected_2 = 0;
n_seg_selected_cnt = 0;
n_seg_selected_cnt_2 = 0;
Nchar = 0;

if phs == 1
    NKR = 2*N_reads_max/DSF;
else
    NKR = N_reads_max/DSF;
end
for kr = 1:1:NKR
    if phs == 1
        % Qmin = -1;
        % while( Qmin < 0 )
            if R_mode == 0
                r_or_c = 0;
                prev_id = kr;
                pei = 0;
                if f_type == 'q'
                    [rd_seq, rd_len, Qa, Qm, Nmin, n_rds] = f04_fastq_get_read_v09( fp_fa, ...
                        cfg.min_read_length, cfg.min_read_length-Min_tail_length, Min_tail_length, ...
                        Normalized_Dist_threshold, Q_ref, DSF, -1, cfg.nominal_read_length, cfg.bool_read_filter );
                    Cdep = 1;
                    % N_reads_sofar = N_reads_sofar + n_rds;
                    % Qmin = round(Qm);
                else
                    [rd_seq, rd_len, Cdep, Qa, Qmin, Nmin, k_drop_time] = ...
                        f04_fasta_get_read_v09( fp_fa, cfg.min_read_length, DSF, 0, Min_tail_length, ...
                        Normalized_Dist_threshold, cfg.nominal_read_length, cfg.bool_read_filter );
                    Cdep = 1;
                    % N_reads_sofar = N_reads_sofar + 1;
                end
            else
                r_or_c = 0;
                prev_id = floor( (kr+1)/2 );
                pei = mod(kr,2);
                
                if f_type == 'q'
                    if pei == 0
                        fpt = fp_fa;
                    else
                        fpt = fp_fp;
                    end
                    [rd_seq, rd_len, Qa, Qm, Nmin, n_rds] = f04_fastq_get_read_v09( fpt, ...
                        cfg.min_read_length, cfg.min_read_length-Min_tail_length, Min_tail_length, ...
                        Normalized_Dist_threshold, Q_ref, DSF, -1, cfg.nominal_read_length, cfg.bool_read_filter );
                    Cdep = 1;
                    % N_reads_sofar = N_reads_sofar + 1;
                else
                    if pei == 0
                        fpt = fp_fa;
                    else
                        fpt = fp_fp;
                    end
                    [rd_seq, rd_len, Cdep, Qa, Qmin, Nmin, k_drop_time] = ...
                        f04_fasta_get_read_v09( fpt, cfg.min_read_length, DSF, 0, Min_tail_length, ...
                        Normalized_Dist_threshold, cfg.nominal_read_length, cfg.bool_read_filter );
                    Cdep = 1;
                    % N_reads_sofar = N_reads_sofar + 1;
                end
            end
            if rd_len < 0 || N_reads_sofar > NKR % N_reads_max/DSF
                % break;
            else
                N_reads_sofar = N_reads_sofar + 1;
            end
        % end
    else
        if Mode_Rd == 0
            [rd_seq, rd_len, Cdep, r_or_c, prev_id, pei, k_drop_time] = f04_fasta_get_read_v09( fp_fa, ...
                cfg.min_read_length, DSF, 1, 0, Normalized_Dist_threshold, cfg.nominal_read_length, 0 );
            N_reads_sofar = kr; %N_reads_sofar + 1;
        else
            r_or_c = Reads_all(kr).r_or_c;
            prev_id = Reads_all(kr).prev_id;
            pei = Reads_all(kr).pei;
            rd_len = Reads_all(kr).len;
            rd_seq = Reads_all(kr).seq(1:rd_len);
            Cdep = round(Reads_all(kr).dep);
            k_drop_time = Reads_all(kr).k_drop;
            N_reads_sofar = kr; %N_reads_sofar + 1;
        end
    end
    
    if rd_len < 0 || N_reads_sofar > NKR % N_reads_max/DSF
        break;
    end
%     if phs == 1 % && R_mode > 0
%         fprintf( fp_s(pei+1), '>%d\n', prev_id );
%         if rd_len > 0
%             fprintf( fp_s(pei+1), '%s\n', sub_NumSeq2NTstr( rd_seq ) );
%         else
%             fprintf( fp_s(pei+1), '\n' );
%         end
%     end
    if rd_len < cfg.min_read_length % || Qa < Q_ave || Qmin < Q_min || (Q_ref > 0 && Nmin > N_min)
        % discard the read
        read_cntg_map( prev_id, pei+1 ) = -1;
    else
        if N_reads_sofar >= nr_time_start && xcntr == 0
            xcntr = 1;
            kclk1 = clock();
        end
        k = k + 1;
        k_disp_flag = 0;
        if phs == 1
            N_reads_valid = N_reads_valid + round(Cdep);
        end
        N_reads_processed = N_reads_processed + round(Cdep);
        
        if rd_len > Min_ovlp_dep
            Ave_num_error_normalized_2 = ceil( rd_len*Normalized_Dist_threshold );
            Seg_select_threshold_2 = round( max( (rd_len - (sim_prm.num_N_mer-Seg_sel_dsr))/Seg_sel_dsr - Ave_num_error_normalized_2*(sim_prm.num_N_mer/Seg_sel_dsr), min_sel_th ) );
        else
            Seg_select_threshold_2 = Seg_select_threshold;
        end
        
        if n_seg_found > 0
            [rd_seq_Nmer1, rd_len_Nmer1] = sub_NumSeq2NmerSeq( rd_seq(1:rd_len), sim_prm.num_N_mer );
            priority(1:n_seg_maxid) = sum( Seg_Nmer_map( 1:n_seg_maxid, rd_seq_Nmer1(1:Seg_sel_dsr:rd_len_Nmer1)+1 ), 2 ); 
        end
        
        sel_seg_idx2 = find( priority(1:n_seg_maxid) >= Seg_select_threshold_2 );
        sel_seg_idx = sel_seg_idx2( b_seg_valid(sel_seg_idx2) > 0 );
        n_seg_considered = length( sel_seg_idx );
        
        dst_all = 10000.*ones(n_seg_considered,1);
        pos_all = zeros(n_seg_considered,1);
        ref_idx = zeros(n_seg_considered,1);
        ovlp_depth = zeros(n_seg_considered,1);
        b_em = zeros(n_seg_considered,1);
        b_po = zeros(n_seg_considered,1);
        bsv = b_seg_valid( sel_seg_idx );
        pri = priority( sel_seg_idx );
        
        if n_seg_considered > 0
            for mm = 1:n_seg_considered
                if bsv(mm) > 0 && pri(mm) >= Seg_select_threshold_2 %b_seg_select(m) > 0
                    % check if read is embedded in the current segment
                    ssi = sel_seg_idx(mm);
                    % [pos_all(mm), ref_idx(mm), dst_all(mm)] = sub_tco( seg(ssi).seq(1:seg(ssi).len), rd_seq, cfg.norm_dist_threshold );
                    [pos_all(mm), ref_idx(mm), dst_all(mm)] = sub_tco( seg(ssi).seq(1:seg(ssi).len), rd_seq, cfg.norm_dist_threshold, cfg.bool_ss_ind );
                    if pos_all(mm) ~= 0
                        % if pos < 0, the 1st seq is embedded in the 2nd
                        % Here, it must be always that pos >= 0
                        b_em(mm) = 1;
                        if dst_all(mm) == 0
                            break;
                        end
                    end
                end
            end % for m = 1:1:M
            n_seg_selected = n_seg_selected + n_seg_considered; %sum(b_seg_select(1:n_seg_maxid));
            n_seg_selected_cnt = n_seg_selected_cnt + 1; %sum(b_seg_select(1:n_seg_maxid));
        end
        if sum( b_em ) == 0
            
            if cfg.b_nmer_map_s > 0 && n_seg_found > 0
                priority(1:n_seg_maxid) = sum( Seg_Nmer_map2( 1:n_seg_maxid, rd_seq_Nmer1(1:Seg_sel_dsr:rd_len_Nmer1)+1 ), 2 ); 
            end
            sel_seg_idx1 = find( priority(1:n_seg_maxid) >= Seg_select_threshold );
            sel_seg_idx = sel_seg_idx1( b_seg_valid(sel_seg_idx1) > 0 );
            n_seg_considered = length( sel_seg_idx );

            dst_all = 10000.*ones(n_seg_considered,1);
            pos_all = zeros(n_seg_considered,1);
            ref_idx = zeros(n_seg_considered,1);
            ovlp_depth = zeros(n_seg_considered,1);
            b_po = zeros(n_seg_considered,1);
            bsv = b_seg_valid( sel_seg_idx );
            pri = priority( sel_seg_idx );
            
            for mm = 1:n_seg_considered
                if bsv(mm) > 0 && pri(mm) >= Seg_select_threshold %b_seg_select(m) > 0
                    % check if read is embedded in the current segment
                    if rd_len > Min_ovlp_dep
                        ssi = sel_seg_idx(mm);
                        % check if current read is partially overlapped with the current segment
                        % [pos_all(mm), ovlp_depth(mm), dst_all(mm)] = sub_tpo( seg(ssi).seq(1:seg(ssi).len), rd_seq, cfg_val );
                        [pos_all(mm), ovlp_depth(mm), dst_all(mm)] = sub_tpo( seg(ssi).seq(1:seg(ssi).len), rd_seq, cfg_val );
                        if pos_all(mm) ~= 0
                            % if abs(pos_all(m))+abs(ovlp_depth(m)) ~= rd_len
                            %     disp(' ERROR: f04_check_partial_overlap_v18' );
                            % end
                            % pos < 0 means rd_seq overlaps with seg_seq at its end
                            % ovlp_dep < 0 means the complement of rd_seq overlaps with seg_seq
                            b_po(mm) = 1;
                        end
                    end
                end
            end % for m = 1:1:M
            n_seg_selected_2 = n_seg_selected_2 + n_seg_considered; %sum(b_seg_select(1:n_seg_maxid));
            n_seg_selected_cnt_2 = n_seg_selected_cnt_2 + 1; %sum(b_seg_select(1:n_seg_maxid));
        end
        
        n_partial_ovlp = sum( b_po );
        n_embedded = sum( b_em );

        % collect stats
        n_em_clipped = min( n_embedded, 3 )+1;
        n_po_clipped = min( n_partial_ovlp, 3 )+1;
        n_em_po_pdf( n_em_clipped, n_po_clipped ) = n_em_po_pdf( n_em_clipped, n_po_clipped ) + 1;

        if n_embedded == 0 

            switch( n_partial_ovlp )

            case 0,
                % This read doesn't match any segments in the library
                % then, add it to the segment libray
                if n_seg_found < Max_pipeline 
                    % select an empty slot
                    n_seg_found = n_seg_found + 1;
                    idxs_tmp = find( b_seg_valid(1:n_seg_maxid) == uint8(0), 1, 'first' );
                    if isempty(idxs_tmp)
                        n_seg_maxid = n_seg_maxid + 1;
                        d = n_seg_maxid;
                    else
                        d = idxs_tmp(1);
                    end
                    b_seg_valid(d) = uint8(1);
                    idx = d;
                    
                    % Initialize segment info
                    if r_or_c == 0
                        cntg_cnt = cntg_cnt + 1;
                        cntg_valid_ind(cntg_cnt) = 0;
                        read_cntg_map( prev_id, pei+1 ) = cntg_cnt;
                        seg(idx).id = cntg_cnt; 
                    else
                        seg(idx).id = prev_id; 
                    end
                    seg(idx).stm = k; %seg_stm(idx) = k;
                    seg(idx).len = rd_len; %seg_len(idx) = rd_len;
                    seg(idx).lin = rd_len; %seg_lin(idx) = rd_len;
                    seg_cvg_tmp2(:,1:rd_len) = 0;
                    seg_cvg_tmp2( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    if seg(idx).blen < rd_len
                        while( seg(idx).blen < rd_len )
                            seg(idx).blen = seg(idx).blen + Initial_buf_length;
                        end
                        seg(idx).seq = zeros(1,seg(idx).blen, 'int8');
                        seg(idx).cvg_dep = zeros(4,seg(idx).blen,'uint32');
                    end
                    seg(idx).seq(1:rd_len) = rd_seq; %seg_seq_lib(idx,1:seg_len(idx)) = rd_seq; 
                    seg(idx).cvg_dep(:,1:rd_len) = seg_cvg_tmp2(:,1:rd_len); %f03_set_cvg( seg_seq_lib(idx,1:seg_len(idx)) );
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                    % End: Initialize segment info

                    [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                    [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                    Seg_Nmer_map( idx, : ) = 0;
                    Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;

                    if cfg.b_nmer_map_s > 0
                        Seg_Nmer_map2( idx, : ) = 0;
                        if seg(idx).len <= cfg.nominal_read_length
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                        else
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                        end
                    end
                    if sim_prm.b_inter_disp > 0
                        str_disp = sprintf('N_seg: %4d/%4d --> New seg [# %2d] of LEN %4d added to library', ...
                            n_seg_considered, n_seg_found, idx, seg(idx).len );
                        disp(str_disp);
                    end
                    
                else
                    if N_dropped == 0;
                        k_first_drop = k;
                    end
                    N_dropped = N_dropped +1;
                    N_reads_processed = N_reads_processed - ( round(Cdep) );
                    if Mode_Rd == 0
                        fprintf(fp_fd,'> %f\t%d\t%d\t%d\t%d\n', Cdep, k, r_or_c, prev_id, pei ); 
                        fprintf(fp_fd, '%s\n', sub_NumSeq2NTstr(rd_seq) );  
                    else
                        Dropped(N_dropped).r_or_c = r_or_c;
                        Dropped(N_dropped).prev_id = prev_id;
                        Dropped(N_dropped).pei = pei;
                        Dropped(N_dropped).len = rd_len;
                        Dropped(N_dropped).seq = rd_seq;
                        Dropped(N_dropped).dep = Cdep;
                        Dropped(N_dropped).k_drop = k;
                    end
                end

                if n_seg_found >= Max_pipeline && Max_pipeline < cfg.max_pipeline
                    % fprintf('Allocating more memories ... ' );
                    Max_contig_length = max( cfg.max_contig_length, max(double([seg(1:n_seg_maxid).len]).*double(b_seg_valid(1:n_seg_maxid))) );
                    Max_seg_length_for_combining = round(Max_seg_length_for_combining*Max_contig_length/cfg.max_contig_length);
                    % b_seg_valid = [b_seg_valid zeros(1,cfg.num_pipeline, 'int16')];
                    % priority = [priority; zeros(cfg.num_pipeline,1, 'int16')];
                    % Max_pipeline = Max_pipeline + cfg.num_pipeline;
                    % seg = [seg; repmat( sg, cfg.num_pipeline, 1 )];
                    % fprintf('done.\n' );
                end

            case 1,
                % n_partial_ovlp > 0
                % This read match two segments in the library
                % connect the two segment and update library

                sel_tmp = b_po.*abs(ovlp_depth);
                [axx, idx1_t] = max( abs( sel_tmp ) );
                idx1 = sel_seg_idx(idx1_t);

                idx = idx1;
                pos = pos_all(idx1_t);
                ovlp_dep = ovlp_depth(idx1_t);
                if pos < 0
                    pos = -pos;
                    rd_seq = 3-rd_seq(end:-1:1);
                end

                % Extend the segment
                % Note: abs(ovlp_depth) + abs(pos) = Seg Len
                if ovlp_dep > 0 % segment front
                    tlen = seg(idx).len + pos;
                    if seg(idx).blen < tlen
                        seg_cvg_tmp2(:,1:seg(idx).len) = seg(idx).cvg_dep(:,1:seg(idx).len);
                        while( seg(idx).blen < tlen )
                            seg(idx).blen = seg(idx).blen + Initial_buf_length;
                        end
                        seg(idx).seq = zeros(1,seg(idx).blen, 'int8');
                        seg(idx).cvg_dep = zeros(4,seg(idx).blen,'uint32');
                        seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp2(:,1:seg(idx).len);
                    end
                    seg(idx).cvg_dep(:,1:tlen) = [ zeros(4,pos,'uint32')  seg(idx).cvg_dep(:,1:seg(idx).len) ];
                    seg_cvg_tmp2(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp2( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg(idx).cvg_dep(:,1:rd_len) = ...
                        seg(idx).cvg_dep(:,1:rd_len) + seg_cvg_tmp2(:,1:rd_len); %f03_set_cvg( rd_seq );
                    seg(idx).len = tlen; %seg(idx).len + pos; %seg_len(idx) = seg_len(idx) + pos;
                    [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                    seg(idx).seq = int8( mxi(1:seg(idx).len)-1 );
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                else % segment end
                    tlen = seg(idx).len + pos;
                    if seg(idx).blen < tlen
                        seg_cvg_tmp2(:,1:seg(idx).len) = seg(idx).cvg_dep(:,1:seg(idx).len);
                        while( seg(idx).blen < tlen )
                            seg(idx).blen = seg(idx).blen + Initial_buf_length;
                        end
                        seg(idx).seq = zeros(1,seg(idx).blen, 'int8');
                        seg(idx).cvg_dep = zeros(4,seg(idx).blen,'uint32');
                        seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp2(:,1:seg(idx).len);
                    end
                    % seg(idx).cvg_dep(:,1:tlen) = [ seg(idx).cvg_dep(:,1:seg(idx).len) zeros(4,pos,'uint32') ];
                    seg(idx).cvg_dep(:,seg(idx).len+1:tlen) = zeros(4,pos,'uint32');
                    seg(idx).len = tlen; %seg(idx).len + pos; %seg_len(idx) = seg_len(idx) + pos;
                    seg_cvg_tmp2(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp2( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg(idx).cvg_dep(:,seg(idx).len-rd_len+1:seg(idx).len) = ...
                        seg(idx).cvg_dep(:,seg(idx).len-rd_len+1:seg(idx).len) + seg_cvg_tmp2(:,1:rd_len); %f03_set_cvg( rd_seq );
                    [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                    seg(idx).seq(1:seg(idx).len) = int8( mxi(1:seg(idx).len)-1 );
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                end

                [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                Seg_Nmer_map( idx, : ) = 0;
                Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;
                seg(idx).stm = k;
                if cfg.b_nmer_map_s > 0
                    Seg_Nmer_map2( idx, : ) = 0;
                    if seg(idx).len <= cfg.nominal_read_length
                        Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                    else
                        Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                    end
                end
                if r_or_c == 0
                    read_cntg_map( prev_id, pei+1 ) = seg(idx).id;
                else
                    cntg_valid_ind(prev_id) = seg(idx).id; 
                end

            otherwise,
                %if n_partial_ovlp > 1
                sel_tmp = b_po.*abs(ovlp_depth);
                [axx, idx1_t] = max( abs( sel_tmp ) );
                sel_tmp(idx1_t) = 0;
                idx1 = sel_seg_idx(idx1_t);

                idx = idx1;
                pos = pos_all(idx1_t);
                ovlp_dep = ovlp_depth(idx1_t);
                if pos < 0
                    pos = -pos;
                    rd_seq = 3-rd_seq(end:-1:1);
                end

                % Extend the segment
                seg_seq_tmp1(1:seg(idx).len) = seg(idx).seq(1:seg(idx).len); 
                seg_cvg_tmp1(:,1:seg(idx).len) = seg(idx).cvg_dep(:,1:seg(idx).len);
                % Note: abs(ovlp_depth) + abs(pos) = Seg Len
                if ovlp_dep > 0 % segment front
                    seg_len_tmp1 = seg(idx).len + pos;
                    seg_cvg_tmp1(:,1:seg(idx).len+pos) = [ zeros(4,pos,'uint32')  seg_cvg_tmp1(:,1:seg(idx).len) ];
                    seg_cvg_tmp2(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp2( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg_cvg_tmp1(:,1:rd_len) = ...
                        seg_cvg_tmp1(:,1:rd_len) + seg_cvg_tmp2(:,1:rd_len); %f03_set_cvg( rd_seq );
                    [mxv(1:seg_len_tmp1), mxi(1:seg_len_tmp1)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1) );
                    seg_seq_tmp1(1:seg_len_tmp1) = int8( mxi(1:seg_len_tmp1)-1 );
                else % segment end
                    seg_len_tmp1 = seg(idx).len+pos;
                    seg_cvg_tmp1(:,1:seg_len_tmp1) = [ seg_cvg_tmp1(:,1:seg(idx).len) zeros(4,pos,'uint32') ];
                    seg_cvg_tmp2(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp2( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg_cvg_tmp1(:,seg_len_tmp1-rd_len+1:seg_len_tmp1) = ...
                        seg_cvg_tmp1(:,seg_len_tmp1-rd_len+1:seg_len_tmp1) + seg_cvg_tmp2(:,1:rd_len); %f03_set_cvg( rd_seq );
                    [mxv(1:seg_len_tmp1), mxi(1:seg_len_tmp1)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1) );
                    seg_seq_tmp1(1:seg_len_tmp1) = int8( mxi(1:seg_len_tmp1)-1 );
                end

                seg(idx).len = seg_len_tmp1;
                tlen = seg(idx).len;
                if seg(idx).blen < tlen
                    while( seg(idx).blen < tlen )
                        seg(idx).blen = seg(idx).blen + Initial_buf_length;
                    end
                    seg(idx).seq = zeros(1,seg(idx).blen, 'int8');
                    seg(idx).cvg_dep = zeros(4,seg(idx).blen,'uint32');
                end
                seg(idx).seq(1:seg(idx).len) = seg_seq_tmp1(1:seg(idx).len); 
                seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp1(:,1:seg(idx).len);
                seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                
                [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                Seg_Nmer_map( idx, : ) = 0;
                Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;
                seg(idx).stm = k;
                if cfg.b_nmer_map_s > 0
                    Seg_Nmer_map2( idx, : ) = 0;
                    if seg(idx).len <= cfg.nominal_read_length
                        Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                    else
                        Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                    end
                end
                if r_or_c == 0
                    read_cntg_map( prev_id, pei+1 ) = seg(idx).id;
                else
                    cntg_valid_ind(prev_id) = seg(idx).id; 
                end
                    
                if n_partial_ovlp > 1
                    for kt = 2:1:n_partial_ovlp
                        [axx, idx2_t] = max( abs( sel_tmp ) );
                        sel_tmp(idx2_t) = 0;
                        idx2 = sel_seg_idx(idx2_t);

                        if max( seg(idx1).len, seg(idx2).len ) < Max_seg_length_for_combining % scfg.b_seg_combining == 0
                        
                            % [pos, ovlp_dep, axx] = sub_tpo( seg_seq_tmp1(:,1:seg_len_tmp1), seg(idx2).seq(1:seg(idx2).len), cfg_val );
                            [pos, ovlp_dep, dst] = sub_tpo( seg_seq_tmp1(:,1:seg_len_tmp1), seg(idx2).seq(1:seg(idx2).len), cfg_val );
                            %pos = 0; % switch
                            if pos ~= 0
                                
                                seg_len_tmp2 = seg(idx2).len;
                                seg_seq_tmp2(1:seg_len_tmp2) = seg(idx2).seq(1:seg(idx2).len);
                                seg_cvg_tmp2(:,1:seg_len_tmp2) = seg(idx2).cvg_dep(:,1:seg(idx2).len);
                                
                                if abs(pos)+abs(ovlp_dep) ~= seg_len_tmp2
                                    str = sprintf(' ERROR: f04_check_partial_overlap_v18 - %d', dst );
                                    disp( str );
                                end
                                if pos < 0
                                    pos = -pos;
                                    seg_seq_tmp2(1:seg_len_tmp2) = 3-seg_seq_tmp2(seg_len_tmp2:-1:1);
                                    seg_cvg_tmp2(:,1:seg_len_tmp2) = seg_cvg_tmp2(4:-1:1,seg_len_tmp2:-1:1);
                                end

                                if ovlp_dep > 0 % segment front
                                    seg_cvg_tmp1(:,1:seg_len_tmp1+pos) = [ zeros(4,pos,'uint32')  seg_cvg_tmp1(:,1:seg_len_tmp1) ];
                                    seg_cvg_tmp1(:,1:seg_len_tmp2) = seg_cvg_tmp1(:,1:seg_len_tmp2) + seg_cvg_tmp2(:,1:seg_len_tmp2);
                                    [mxv(1:seg_len_tmp1+pos), mxi(1:seg_len_tmp1+pos)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1+pos) );
                                    seg_seq_tmp1(1:seg_len_tmp1+pos) = int8( mxi(1:seg_len_tmp1+pos)-1 );
                                    seg_len_tmp1 = seg_len_tmp1+pos;
                                else % segment end
                                    seg_cvg_tmp1(:,1:seg_len_tmp1+pos) = [ seg_cvg_tmp1(:,1:seg_len_tmp1) zeros(4,pos,'uint32') ];
                                    seg_len_tmp1 = seg_len_tmp1+pos;
                                    seg_cvg_tmp1(:,seg_len_tmp1-seg_len_tmp2+1:seg_len_tmp1) = ...
                                        seg_cvg_tmp1(:,seg_len_tmp1-seg_len_tmp2+1:seg_len_tmp1) + seg_cvg_tmp2(:,1:seg_len_tmp2);
                                    [mxv(1:seg_len_tmp1), mxi(1:seg_len_tmp1)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1) );
                                    seg_seq_tmp1(1:seg_len_tmp1) = int8( mxi(1:seg_len_tmp1)-1 );
                                end

                                idx = min( idx1, idx2 );
                                seg(idx).len = seg_len_tmp1;
                                tlen = seg(idx).len;
                                if seg(idx).blen < tlen
                                    while( seg(idx).blen < tlen )
                                        seg(idx).blen = seg(idx).blen + Initial_buf_length;
                                    end
                                    seg(idx).seq = zeros(1,seg(idx).blen, 'int8');
                                    seg(idx).cvg_dep = zeros(4,seg(idx).blen,'uint32');
                                end
                                seg(idx).seq(1:seg(idx).len) = seg_seq_tmp1(1:seg(idx).len); 
                                seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp1(:,1:seg(idx).len);
                                seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );

                                [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                                [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                                Seg_Nmer_map( idx, : ) = 0;
                                Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;
                                seg(idx).stm = k;
                                if cfg.b_nmer_map_s > 0
                                    Seg_Nmer_map2( idx, : ) = 0;
                                    if seg(idx).len <= cfg.nominal_read_length
                                        Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                                    else
                                        Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                                    end
                                end
                                idx_tmp = idx;

                                idx = max( idx1, idx2 );
                                cntg_valid_ind( seg(idx).id ) = seg(idx_tmp).id; 
                                
                                b_seg_valid(idx) = 0;
                                seg(idx) = sg;
                                n_seg_found = n_seg_found - 1;
                                
                                break;
                            end
                        else
                        end % cfg.b_seg_combining == 0
                    end % for kt
                end % seg_len(idx1) < cfg.max_seg_length_for_combining
                
                if sim_prm.b_inter_disp > 0
                    str_disp = sprintf('Nseg:%4d%4d, Seg#:%4d, Len:%5d, k:%6d/%6d, N_se:(%2d, %2d)', ...
                        n_seg_considered, n_seg_found, idx1, seg_len_tmp1, k, N_reads_max/DSF, ...
                        sum(n_em_po_pdf(4,:)), sum(n_em_po_pdf(:,4)) );
                    disp(str_disp);
                end

            end % switch

        else % n_embedded > 0 

            if n_embedded == 1

                [mxvx, idx_t] = max( b_em );
                pos = pos_all(idx_t);
                % dist = dst_all(idx_t);
                ref = ref_idx(idx_t);
                idx = sel_seg_idx(idx_t);
                if pos < 0
                    pos = -pos;
                    if ref == 0
                        rd_seq = 3-rd_seq(end:-1:1);
                    else
                        seg(idx).cvg_dep(:,1:seg(idx).len) = seg(idx).cvg_dep(4:-1:1,seg(idx).len:-1:1);
                        seg(idx).seq(1:seg(idx).len) = 3-seg(idx).seq(seg(idx).len:-1:1);
                    end
                end

                % Initializing variables
                % Update seg coverage and seq, Mark read as aligned
                if ref == 0
                    seg_cvg_tmp1(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp1( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg(idx).cvg_dep(:,pos:1:pos+rd_len-1) = ...
                        seg(idx).cvg_dep(:,pos:1:pos+rd_len-1) + seg_cvg_tmp1(:,1:rd_len); %f03_set_cvg( rd_seq_tmp(1:rd_len_tmp) );
                    % if dist > 0
                        [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                        seg(idx).seq(1:seg(idx).len) = int8( mxi(1:seg(idx).len)-1 );
                        seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                    % end
                else
                    seg_cvg_tmp1(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp1( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg_cvg_tmp1(:,pos:1:pos+seg(idx).len-1) =  seg_cvg_tmp1(:,pos:1:pos+seg(idx).len-1) +...
                        seg(idx).cvg_dep(:,1:seg(idx).len);
                    seg(idx).len = rd_len;
                    seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp1(:,1:seg(idx).len); %seg_cvg_tmp2(:,1:seg(idx).len);
                    [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                    seg(idx).seq(1:seg(idx).len) = int8( mxi(1:seg(idx).len)-1 );
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );

                    [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                    [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                    Seg_Nmer_map( idx, : ) = 0;
                    Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;
                    % seg(idx).stm = k;
                    if cfg.b_nmer_map_s > 0
                        Seg_Nmer_map2( idx, : ) = 0;
                        if seg(idx).len <= cfg.nominal_read_length
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                        else
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                        end
                    end
                end
                % End: Update seg coverage and seq, Mark read as aligned
                if r_or_c == 0
                    read_cntg_map( prev_id, pei+1 ) = seg(idx).id;
                else
                    cntg_valid_ind(prev_id) = seg(idx).id; 
                end

                % display status
                if sim_prm.b_inter_disp > 0
                    Nlen = seg(idx).len;
                    str_disp = sprintf('Nseg:%4d/%4d, Seg#:%4d, Len:%5d, k:%6d/%6d, N_se:(%2d, %2d)', ...
                        n_seg_considered, n_seg_found, idx, Nlen, k, N_reads_max/DSF, ...
                        sum(n_em_po_pdf(4,:)), sum(n_em_po_pdf(:,4)) );
                    disp( str_disp );
                end
                
            else % n_embedded > 1

                sel_tmp = double( dst_all.*b_em + (1-b_em).*10000 ); 
                [mnvx, idx1_t] = min( sel_tmp );
                sel_tmp(idx1_t) = 10000;
                [mxvx, idx2_t] = min( sel_tmp );
                idx1 = sel_seg_idx(idx1_t);
                idx2 = sel_seg_idx(idx2_t);
                
                if seg(idx1).len > seg(idx2).len
                    idx_l = idx1;
                    idx_s = idx2;
                    idx_l_t = idx1_t;
                else
                    idx_l = idx2;
                    idx_s = idx1;
                    idx_l_t = idx2_t;
                end
                clear idx1;
                clear idx2;

                idx = idx_l;
                idx_t = idx_l_t;
                pos = pos_all(idx_t);
                % dist = dst_all(idx_t);
                ref = ref_idx(idx_t);
                if pos < 0
                    pos = -pos;
                    if ref == 0
                        rd_seq = 3-rd_seq(end:-1:1);
                    else
                        seg(idx).cvg_dep(:,1:seg(idx).len) = seg(idx).cvg_dep(4:-1:1,seg(idx).len:-1:1);
                        seg(idx).seq(1:seg(idx).len) = 3-seg(idx).seq(seg(idx).len:-1:1);
                    end
                end

                % Initializing variables
                % Update seg coverage and seq, Mark read as aligned
                if ref == 0
                    seg_cvg_tmp1(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp1( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg(idx).cvg_dep(:,pos:1:pos+rd_len-1) = ...
                        seg(idx).cvg_dep(:,pos:1:pos+rd_len-1) + seg_cvg_tmp1(:,1:rd_len); %f03_set_cvg( rd_seq_tmp(1:rd_len_tmp) );
                    % if dist > 0
                        [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                        seg(idx).seq(1:seg(idx).len) = int8( mxi(1:seg(idx).len)-1 );
                        seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );
                    % end
                else
                    seg_cvg_tmp1(:,1:rd_len) = uint32(0);
                    seg_cvg_tmp1( (0:4:4*(rd_len-1))+double(rd_seq)+1 ) = round(Cdep);
                    seg_cvg_tmp1(:,pos:1:pos+seg(idx).len-1) =  seg_cvg_tmp1(:,pos:1:pos+seg(idx).len-1) +...
                        seg(idx).cvg_dep(:,1:seg(idx).len);
                    seg(idx).len = rd_len;
                    seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp1(:,1:seg(idx).len); 
                    [mxv(1:seg(idx).len), mxi(1:seg(idx).len)] = max( seg(idx).cvg_dep(:,1:seg(idx).len) );
                    seg(idx).seq(1:seg(idx).len) = int8( mxi(1:seg(idx).len)-1 );
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );

                    [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), sim_prm.num_N_mer );
                    [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), sim_prm.num_N_mer );
                    Seg_Nmer_map( idx, : ) = 0;
                    Seg_Nmer_map( idx, [seg_seq_Nmer1(1:sl1); seg_seq_Nmer2(1:sl2)]+1 ) = 1;
                    % seg(idx).stm = k;
                    if cfg.b_nmer_map_s > 0
                        Seg_Nmer_map2( idx, : ) = 0;
                        if seg(idx).len <= cfg.nominal_read_length
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1 seg_seq_Nmer2]+1 ) = 1;
                        else
                            Seg_Nmer_map2( idx, [seg_seq_Nmer1(1:slen) seg_seq_Nmer1(end-slen+1:end) seg_seq_Nmer2(1:slen) seg_seq_Nmer2(end-slen+1:end)]+1 ) = 1;
                        end
                    end
                end
                if r_or_c == 0
                    read_cntg_map( prev_id, pei+1 ) = seg(idx).id;
                else
                    cntg_valid_ind(prev_id) = seg(idx).id; 
                end
                % End: Update seg coverage and seq

                if seg(idx_s).len > 0
                    [pos, ref, dst] = sub_tco( seg(idx_l).seq(1:seg(idx_l).len), seg(idx_s).seq(1:seg(idx_s).len), cfg.norm_dist_threshold, cfg.bool_ss_ind );
                else
                    pos = 0;
                end
                if pos == 0
                    if sim_prm.b_inter_disp > 0
                        str = sprintf('Two seg NOT merged ... %d', dst );
                        disp(str);
                    end
                else
                    idx = idx_l;
                    seg_len_tmp1 = seg(idx).len;
                    seg_seq_tmp1(1:seg(idx).len) = seg(idx).seq(1:seg(idx).len); 
                    seg_cvg_tmp1(:,1:seg(idx).len) = seg(idx).cvg_dep(:,1:seg(idx).len);

                    idx = idx_s;
                    seg_len_tmp2 = seg(idx).len;
                    seg_seq_tmp2(1:seg(idx).len) = seg(idx).seq(1:seg(idx).len); 
                    seg_cvg_tmp2(:,1:seg(idx).len) = seg(idx).cvg_dep(:,1:seg(idx).len);
                
                    if pos < 0
                        pos = -pos;
                        if ref == 0
                            seg_cvg_tmp2(:,1:seg_len_tmp2) = seg_cvg_tmp2(4:-1:1,seg_len_tmp2:-1:1);
                            seg_seq_tmp2(1:seg_len_tmp2) = 3-seg_seq_tmp2(seg_len_tmp2:-1:1);
                        else
                            seg_cvg_tmp1(:,1:seg_len_tmp1) = seg_cvg_tmp1(4:-1:1,seg_len_tmp1:-1:1);
                            seg_seq_tmp1(1:seg_len_tmp1) = 3-seg_seq_tmp1(seg_len_tmp1:-1:1);
                        end
                    end
                    if ref == 0
                        seg_cvg_tmp1(:,pos:1:pos+seg_len_tmp2-1) = ...
                            seg_cvg_tmp1(:,pos:1:pos+seg_len_tmp2-1) + seg_cvg_tmp2(:,1:seg_len_tmp2);
                        [mxv(1:seg_len_tmp1), mxi(1:seg_len_tmp1)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1) );
                        seg_seq_tmp1(1:seg_len_tmp1) = int8( mxi(1:seg_len_tmp1)-1 );
                    else
                        seg_cvg_tmp2(:,pos:1:pos+seg_len_tmp1-1) = ...
                            seg_cvg_tmp2(:,pos:1:pos+seg_len_tmp1-1) + seg_cvg_tmp1(:,1:seg_len_tmp1);
                        seg_seq_tmp2(1:seg_len_tmp2) = f03_seq_est( seg_cvg_tmp2(:,1:seg_len_tmp2) );
                        
                        seg_len_tmp1 = seg_len_tmp2;
                        seg_cvg_tmp1(:,1:seg_len_tmp1) = seg_cvg_tmp2(:,1:seg_len_tmp1);
                        [mxv(1:seg_len_tmp1), mxi(1:seg_len_tmp1)] = max( seg_cvg_tmp1(:,1:seg_len_tmp1) );
                        seg_seq_tmp1(1:seg_len_tmp1) = int8( mxi(1:seg_len_tmp1)-1 );
                    end
                    
                    idx = min( idx_l, idx_s );
                    seg(idx).len = seg_len_tmp1;
                    seg(idx).cvg_dep(:,1:seg(idx).len) = seg_cvg_tmp1(:,1:seg(idx).len);
                    seg(idx).seq(1:seg(idx).len) = seg_seq_tmp1(1:seg(idx).len); 
                    seg(idx).ave_cvg_dep = mean( sum( seg(idx).cvg_dep(:,1:seg(idx).len), 1, 'native' ) );

                    if idx ~= idx_l
                        Seg_Nmer_map( idx_s, : ) = Seg_Nmer_map( idx_l, : );
                        % seg(idx_s).stm = k;
                        if cfg.b_nmer_map_s > 0
                            Seg_Nmer_map2( idx_s, : ) = Seg_Nmer_map2( idx_l, : );
                        end
                    end
                    idx_tmp = idx;

                    idx = max( idx_l, idx_s );
                    cntg_valid_ind( seg(idx).id ) = seg(idx_tmp).id; 
                    
                    b_seg_valid(idx) = 0;
                    seg(idx) = sg;
                    n_seg_found = n_seg_found - 1;
                end
                
            end
        end
    end % if rd_mapped_seq(k,1) == 0

    if mod( k, sim_prm.update_period ) == 0 && k > 0 
        % Get stats and discard segs
        seg_cnt = 0; 
        seg_cnt_l = 0;
        seg_cnt_n = 0;
        seg_cnt_g = 0;
        for n = 1:n_seg_maxid
            if b_seg_valid(n) == 1
                ave_cvg_dep_tmp = seg(n).ave_cvg_dep;
                if cfg.contig_ttl > 0 && seg(n).stm < (k - C_TTL) 
                    if  (ave_cvg_dep_tmp < 1.001) || (seg(n).len == seg(n).lin && ave_cvg_dep_tmp < Min_CV_to_remove) || ...
                            (seg(n).len <= cfg.nominal_read_length && ave_cvg_dep_tmp < Min_CV_to_remove)
                        if k > C_TTL && n_seg_found > cfg.max_pipeline-sim_prm.update_period && Mode == 0
                            b_seg_valid(n) = 0;
                            n_seg_found = n_seg_found - 1;
                            % if phs < N_phs_loop
                                if N_dropped == 0;
                                    k_first_drop = k;
                                end
                                N_dropped = N_dropped +1;
                                N_reads_processed = N_reads_processed - round(ave_cvg_dep_tmp);
                                if Mode_Rd == 0
                                    % fprintf(fp_fd,'> %f\t%d\n', (ave_cvg_dep_tmp), k ); 
                                    fprintf(fp_fd,'> %f\t%d\t%d\t%d\t%d\n', (ave_cvg_dep_tmp), k, 1, seg(n).id, 0 ); 
                                    fprintf(fp_fd, '%s\n', sub_NumSeq2NTstr(seg(n).seq(1:seg(n).len)) );  
                                else
                                    Dropped(N_dropped).r_or_c = 1;
                                    Dropped(N_dropped).prev_id = seg(n).id;
                                    Dropped(N_dropped).pei = 0;
                                    Dropped(N_dropped).len = seg(n).len;
                                    Dropped(N_dropped).seq = seg(n).seq(1:seg(n).len);
                                    Dropped(N_dropped).dep = seg(n).ave_cvg_dep;
                                    Dropped(N_dropped).k_drop = k;
                                end
                            % end
                            seg(n) = sg;
                        end
                    else
                        seg_cnt = seg_cnt + 1;
                        if ave_cvg_dep_tmp <= 1.001 
                            seg_cnt_n = seg_cnt_n + 1;
                        end
                        if seg(n).len == seg(n).lin 
                            seg_cnt_g = seg_cnt_g + 1;
                        end
                        if seg(n).len < Min_ovlp_dep 
                            seg_cnt_l = seg_cnt_l + 1;
                        end
                    end
                else
                    seg_cnt = seg_cnt + 1;
                    if ave_cvg_dep_tmp <= 1.001 
                        seg_cnt_n = seg_cnt_n + 1;
                    end
                    if seg(n).len == seg(n).lin 
                        seg_cnt_g = seg_cnt_g + 1;
                    end
                    if seg(n).len < Min_ovlp_dep 
                        seg_cnt_l = seg_cnt_l + 1;
                    end
                end
            end
        end
        pr_n = round( 100*seg_cnt_n/seg_cnt ); % 세그먼트 커버리지 깊이가 1보다 작거나 같은 세그먼트의 비율 
        pr_g = round( 100*seg_cnt_g/seg_cnt ); % 성장하지 않은 세그먼트의 비율 
    end
        
    sav_period = 50000;
    if mod( k, sav_period ) == 0 && cfg.contig_ats <= 0  
        if phs > 1 && k_drop_time > 0
            seg_cnt = 0;
            max_cntg_len = 0;
            for k1 = 1:1:n_seg_maxid
                b_tmp = (b_seg_valid(k1) > 0) && (seg(k1).stm <= k_drop_time);
                if b_tmp > 0 && seg(k1).len > 0
                    if max_cntg_len < seg(k1).len
                        max_cntg_len = seg(k1).len;
                    end
                    if seg(k1).len > seg(k1).lin || seg(k1).ave_cvg_dep > 1 %= Min_CV_to_remove
                        seg_cnt = seg_cnt +1 ;
                        % fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, N_seg + seg_cnt, N_reads_valid, N_reads_total );
                        fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, seg(k1).id, N_reads_valid, N_reads_total );
                        fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg(k1).seq(1:seg(k1).len) ) );
                        for m2 = 1:1:4
                            for m1 = 1:1:seg(k1).len-1
                                fprintf(fp_t,'%d\t', seg(k1).cvg_dep(m2,m1) );
                            end
                            m1 = seg(k1).len;
                            fprintf(fp_t,'%d\n', seg(k1).cvg_dep(m2,m1) );
                        end
                        n_base = n_base + seg(k1).len;
                        cvg_sum = sum( seg(k1).cvg_dep(:,1:seg(k1).len) );
                        cvg_max = max( seg(k1).cvg_dep(:,1:seg(k1).len) );
                        sum_sum = sum_sum + sum(cvg_sum);
                        sum_max = sum_max + sum(cvg_max);
                        
                        b_seg_valid(k1) = 0;
                        seg(k1) = sg;
                    else % Not grown && seg(k1).ave_cvg_dep < Min_CV_to_remove
                        ave_cvg_dep_tmp = seg(k1).ave_cvg_dep;
                        b_seg_valid(k1) = 0;
                        n_seg_found = n_seg_found - 1;
                        % N_dropped = N_dropped +1;
                        N_reads_processed = N_reads_processed - round(ave_cvg_dep_tmp);
                        N_discarded = N_discarded + round(ave_cvg_dep_tmp);
                        
                        cntg_valid_ind( seg(k1).id ) = -2;
                    end
                end
            end
            N_seg = N_seg + seg_cnt;
            n_seg_found = n_seg_found - seg_cnt;
            if max_cntg_len < 1.1*cfg.nominal_read_length && Mode == 0
                % Mode = 1;
            end
        end
    end
    if mod( k, sav_period ) == 0 && cfg.contig_ats > 0
        if k > 0
            sav_length = cfg.contig_ats;
            seg_cnt = 0;
            for n = 1:n_seg_maxid
                if b_seg_valid(n) == 1 && seg(n).stm < (k - sav_length) 
                    % if  seg(n).len > seg(n).lin && ave_cvg_dep_tmp > Min_CV_to_remove
                        seg_cnt = seg_cnt +1 ;
                        % fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, N_seg + seg_cnt, N_reads_valid, N_reads_total );
                        fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, seg(n).id, N_reads_valid, N_reads_total );
                        fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg(n).seq(1:seg(n).len) ) );
                        for m2 = 1:1:4
                            for m1 = 1:1:seg(n).len-1
                                fprintf(fp_t,'%d\t', seg(n).cvg_dep(m2,m1) );
                            end
                            m1 = seg(n).len;
                            fprintf(fp_t,'%d\n', seg(n).cvg_dep(m2,m1) );
                        end
                        k1 = n;
                        n_base = n_base + seg(k1).len;
                        cvg_sum = sum( seg(k1).cvg_dep(:,1:seg(k1).len) );
                        cvg_max = max( seg(k1).cvg_dep(:,1:seg(k1).len) );
                        sum_sum = sum_sum + sum(cvg_sum);
                        sum_max = sum_max + sum(cvg_max);

                        b_seg_valid(n) = 0;
                        seg(n) = sg;
                    % end
                end
            end
            N_seg = N_seg + seg_cnt;
            n_seg_found = n_seg_found - seg_cnt;
        end
    end
    
    if mod( k, sim_prm.disp_period ) == 0 && k > 0
        if k_disp_flag == 0
            k_disp_flag = 1;
            if N_reads_sofar <= nr_time_start 
               tstr = ' '; 
            else
                kclk2 = clock();
                ktime = round( clk2sec( kclk2 - kclk1 ) );
                rtime_sec = max(N_reads_max/DSF - N_reads_sofar,0)*(ktime/(N_reads_sofar - nr_time_start));
                clk = sec2clk( round(rtime_sec) );
                tstr = sprintf('.. Remaining %dd, %dh, %dm', clk(3), clk(4), clk(5) );
            end
            max_cntg_len = max( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] );
            min_cntg_len = min( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] );
            med_cntg_len = round( median( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] ) );
            str_disp = sprintf('   L%d (%d): Prog:%4.1f%%(%dK), Sel:%3.1f/%3.1f, DR:%4.1f%%(%dK), Fnd:%4d(%d), CL:%4d/%d/%d, Stat:%d/%d/%dK %s', ...
                phs, Min_ovlp_dep, (100*N_reads_sofar/(N_reads_max/DSF)), round(N_reads_sofar/1000), ...
                n_seg_selected/n_seg_selected_cnt, n_seg_selected_2/n_seg_selected_cnt_2, ...
                (100*N_dropped/kr), round(N_dropped/1000), ...
                n_seg_found, N_seg, ...
                max_cntg_len, med_cntg_len, min_cntg_len, ...
                pr_n, pr_g, round(N_discarded/1000), tstr );
            fprintf(repmat('\b', 1, Nchar));
            Nchar = fprintf( '%s', str_disp ); 
        end
    end
end % for k = 1:1:N_reads_max

max_cntg_len = max( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] );
min_cntg_len = min( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] );
med_cntg_len = round( median( [seg( b_seg_valid(1:n_seg_maxid)>0 ).len] ) );
Drop_rate = N_dropped/kr;
str_disp = sprintf('   L%d (%d): Prog:%4.1f%%(%dK), Sel:%3.1f/%3.1f, DR:%4.1f%%(%dK), Fnd:%4d(%d), CL:%4d/%d/%d, Stat:%d/%d/%dK ', ...
    phs, Min_ovlp_dep, (100*N_reads_sofar/(N_reads_max/DSF)), round(N_reads_sofar/1000), ...
    n_seg_selected/n_seg_selected_cnt, n_seg_selected_2/n_seg_selected_cnt_2, ...
    (100*Drop_rate), round(N_dropped/1000), ...
    n_seg_found, N_seg, ...
    max_cntg_len, med_cntg_len, min_cntg_len, ...
    pr_n, pr_g, round(N_discarded/1000) );
fprintf(repmat('\b', 1, Nchar));
fprintf( '%s\n', str_disp ); 
fprintf( fp_log, '%s', str_disp ); 

if Mode_Rd == 0
    fclose(fp_fa);
    if phs == 1
        if R_mode > 0
            fclose(fp_fp);
        end
    else
        pause(2);
        delete( fname_fa );
    end
end

str_disp = sprintf('   Loop %d completed %s', phs, datestr(now) );
fprintf('%s,', str_disp);
fprintf( fp_log, '\n%s', str_disp ); 
if sum_sum > 0
    estimated_error_rate = (sum_sum - sum_max)/sum_sum;
    str_disp = sprintf(' Estimated RER: %f%% ', estimated_error_rate*100 );
else
    str_disp = ' ';
end
fprintf('%s', str_disp);
fprintf( fp_log, '%s', str_disp ); 

%% Save to MAT file
cntg_lst_tmp.len = 0;
cntg_lst_tmp.seq = sub_NumSeq2NTstr( [] );
cntg_lst_tmp.consensus = uint16( zeros(4,0) );

if phs == 1
    N_reads_total = N_reads_sofar;
%     fclose( fp_s(1) );
%     if R_mode > 0
%         fclose( fp_s(2) );
%     end
end

% Check if any segs has grown
n_tmp = 0;
for k1 = 1:1:n_seg_maxid
    if b_seg_valid(k1) > 0 && seg(k1).len > seg(k1).lin
        n_tmp = n_tmp + 1;
    else
    end
end
% Stop loop if no segs has grown
if N_dropped == 0 || phs == N_phs_loop || n_tmp < 20 
    
    n_seg_found = n_seg_found - seg_cnt;
    str_disp = sprintf(' .. (%d, %d, %d)', N_seg, n_seg_found, N_dropped );
    fprintf( str_disp );
    fprintf( fp_log, '%s ', str_disp ); 

    str_disp = sprintf('\n   Saving ' );
    fprintf( str_disp );
    fprintf( fp_log, '%s', str_disp ); 

    kmod = round(n_seg_maxid/6);
    N_reads_valid = N_reads_valid - N_dropped - N_discarded;
    seg_cnt = 0;
    for k1 = 1:1:n_seg_maxid
        if b_seg_valid(k1) > 0 
            if seg(k1).ave_cvg_dep >= 1 && seg(k1).len >= seg(k1).lin % cfg.nominal_read_length
                seg_cnt = seg_cnt +1 ;
                % fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, N_seg + seg_cnt, N_reads_valid, N_reads_total );
                fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, seg(k1).id, N_reads_valid, N_reads_total );
                fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg(k1).seq(1:seg(k1).len) ) );
                for m2 = 1:1:4
                    for m1 = 1:1:seg(k1).len-1
                        fprintf(fp_t,'%d\t', seg(k1).cvg_dep(m2,m1) );
                    end
                    m1 = seg(k1).len;
                    fprintf(fp_t,'%d\n', seg(k1).cvg_dep(m2,m1) );
                end
                n_base = n_base + seg(k1).len;
                cvg_sum = sum( seg(k1).cvg_dep(:,1:seg(k1).len) );
                cvg_max = max( seg(k1).cvg_dep(:,1:seg(k1).len) );
                sum_sum = sum_sum + sum(cvg_sum);
                sum_max = sum_max + sum(cvg_max);

                b_seg_valid(k1) = 0;
                seg(k1) = sg;
            else
                n_seg_found = n_seg_found - 1;
                N_dropped = N_dropped +1;
                N_reads_processed = N_reads_processed - round(ave_cvg_dep_tmp);
                if Mode_Rd == 0
                    % fprintf(fp_fd,'> %d\t%d\n', round(ave_cvg_dep_tmp), k ); 
                    fprintf(fp_fd,'> %f\t%d\t%d\t%d\t%d\n', (ave_cvg_dep_tmp), k, 1, seg(k1).id, 0 ); 
                    fprintf(fp_fd, '%s\n', sub_NumSeq2NTstr(seg(k1).seq(1:seg(k1).len)) );  
                else
                    Dropped(N_dropped).r_or_c = 1;
                    Dropped(N_dropped).prev_id = seg(k1).id;
                    Dropped(N_dropped).pei = 0;
                    Dropped(N_dropped).len = seg(k1).len;
                    Dropped(N_dropped).seq = seg(k1).seq(1:seg(k1).len);
                    Dropped(N_dropped).dep = seg(k1).ave_cvg_dep;
                    Dropped(N_dropped).k_drop = k;
                end
                
                b_seg_valid(k1) = 0;
                seg(k1) = sg;
            end
        end
        if mod(k1,kmod) == 0
            fprintf('.');
        end
    end
    N_seg = N_seg + seg_cnt;
    % n_seg_found = n_seg_found - seg_cnt;
    
    if Mode_Rd == 0
        fclose(fp_fd);
    end
    break;
else
    kmod = round(n_seg_maxid/6);
    seg_cnt = 0;
    max_len = 0;
    
    for k1 = 1:1:n_seg_maxid
        if b_seg_valid(k1) > 0 && max_len < seg(k1).len
            max_len = seg(k1).len;
        end
        if phs == 1
            b_tmp = (b_seg_valid(k1) > 0) && (seg(k1).stm <= k_first_drop);
        else
            b_tmp = (b_seg_valid(k1) > 0) && (seg(k1).stm <= k_drop_time);
        end
        if b_tmp > 0 
            if seg(k1).len > seg(k1).lin || seg(k1).ave_cvg_dep >= Min_CV_to_remove
                seg_cnt = seg_cnt +1 ;
                % fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, N_seg + seg_cnt, N_reads_valid, N_reads_total );
                fprintf(fp_t,'>Contig\t%d\t%d\t%d\t%d\n', phs, seg(k1).id, N_reads_valid, N_reads_total );
                fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg(k1).seq(1:seg(k1).len) ) );
                for m2 = 1:1:4
                    for m1 = 1:1:seg(k1).len-1
                        fprintf(fp_t,'%d\t', seg(k1).cvg_dep(m2,m1) );
                    end
                    m1 = seg(k1).len;
                    fprintf(fp_t,'%d\n', seg(k1).cvg_dep(m2,m1) );
                end
                n_base = n_base + seg(k1).len;
                cvg_sum = sum( seg(k1).cvg_dep(:,1:seg(k1).len) );
                cvg_max = max( seg(k1).cvg_dep(:,1:seg(k1).len) );
                sum_sum = sum_sum + sum(cvg_sum);
                sum_max = sum_max + sum(cvg_max);

                b_seg_valid(k1) = 0;
                seg(k1) = sg;
            else 
                
            end
        else
        end
        if mod(k1,kmod) == 0
            fprintf('.');
        end
    end
    N_seg = N_seg + seg_cnt;
    n_seg_found = n_seg_found - seg_cnt;
    str_disp = sprintf(' (%d, %d, %d)', N_seg, n_seg_found, N_dropped );
    disp( str_disp );
    fprintf( fp_log, '%s\n', str_disp ); 
    
    if Mode_Rd == 0
        fclose(fp_fd);
        pause(2);
    end
    if max_len <= cfg.nominal_read_length
        break;
    end
end

if phs > cfg.num_phs_loop
    break;
end
end

if N_dropped < 1
    if Mode_Rd == 0
        delete( fname_fd );
    end
else
    if Mode_Rd == 0
        kmod = round(N_dropped*2/6);
        fname_fa2 = sprintf('%s_tmp%d.fasta', fname_out, phs );
        fp_fa2 = fopen( fname_fa2, 'r' );
        if fp_fa2 < 0
            fprintf('\nFile not exist or cannot open: %s \n', fname_fa );
            return;
        end
        fname_fd = sprintf('%s_dropped.fasta', fname_out );
        fp_fd = fopen( fname_fd, 'w' );
        if fp_fd < 0
            fprintf('\nCannot open: %s \n', fname_fd );
            return;
        end
        n_cnt = 0;
        while(1)
            n_cnt = n_cnt +1;
            line_a = fgets(fp_fa2);
            line_b = fgets(fp_fa2);
            if line_a < 0
                break;
            else
                fprintf( fp_fd, '%s', line_a );
                fprintf( fp_fd, '%s', line_b );
                val = sscanf( line_a(2:end), '%f' );  
                cntg_valid_ind( val(4) ) = -2;
            end
            if mod(n_cnt,kmod) == 0
                fprintf('.');
            end
        end
        fclose(fp_fa2);
        pause(2);
        fclose(fp_fd);
        delete( fname_fa2 );
    else
        kmod = round(N_dropped*2/6);
        fname_fd = sprintf('%s_dropped.fasta', fname_out );
        fp_fd = fopen( fname_fd, 'wt' );
        if fp_fd < 0
            fprintf('\nCannot open: %s \n', fname_fd );
            return;
        end
        for m = 1:1:N_dropped
            fprintf(fp_fd,'> %d\t%d\t%d\t%d\t%d\n', round(Dropped(m).dep), round(Dropped(m).k_drop), Dropped(m).r_or_c, Dropped(m).prev_id, Dropped(m).pei ); 
            fprintf(fp_fd, '%s\n', sub_NumSeq2NTstr(Dropped(m).seq) );  
            if mod(m,kmod) == 0
                fprintf('.');
            end
            cntg_valid_ind( Dropped(m).prev_id ) = -2;
        end
        fclose(fp_fd);
    end
end
fclose(fp_t);

fname_txt = sprintf('%s.rcmap', fname_out );
fp_t = fopen( fname_txt, 'wt' );
fprintf(fp_t, '%d\n', cntg_cnt);
kmod = round(cntg_cnt/6);
for k = 1:1:cntg_cnt
    fprintf(fp_t, '%d\n', cntg_valid_ind(k));
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
Nrt = round( N_reads_total/(R_mode+1) );
fprintf(fp_t, '%d\n', Nrt);
kmod = round(Nrt/6);
for k = 1:1:Nrt
    fprintf(fp_t, '%d\t%d\n', read_cntg_map(k,1), read_cntg_map(k,2) );
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
fclose(fp_t);

str_disp = sprintf(' done. %s (%d, %d, %d)', fname_txt, N_seg, cntg_cnt, max( max( read_cntg_map ) ) );
disp( str_disp );
fprintf( fp_log, '%s\n', str_disp ); 
    
%fprintf('  # Stats.\n' );
fprintf('   Num. of reads total: %d\n', N_reads_total );
fprintf('   Num. of reads processed: %d\n', N_reads_valid );
fprintf('   Num. of reads discarded: %d\n', N_reads_total-N_reads_valid );
fprintf('   Num. of contigs found: %d\n', N_seg );
estimated_error_rate = (sum_sum - sum_max)/sum_sum;
fprintf('   Estimated Read Error Rate: %f%% \n', estimated_error_rate*100 );
fprintf('Contig growing started %s', dstr );
fprintf(' and completed %s\n', datestr(now) );

fprintf(fp_log, '   # Stats.\n' );
fprintf(fp_log, '   Num. of reads total: %d\n', N_reads_total );
fprintf(fp_log, '   Num. of reads processed: %d\n', N_reads_valid );
fprintf(fp_log, '   Num. of reads discarded: %d\n', N_reads_total-N_reads_valid );
fprintf(fp_log, '   Num. of contigs found: %d\n', N_seg );
fprintf(fp_log, '   Estimated Read Error Rate: %f%% \n', estimated_error_rate*100 );
fprintf(fp_log, 'Contig growing started %s', dstr );
fprintf(fp_log, '  and completed %s\n', datestr(now) );
fclose(fp_log);
pause(2);

end % sub_build1_seg_lib_v31_fasta

%% %%%%%%%%% %%
%% functions %%
%% %%%%%%%%% %%

function seq_int8 = f03_seq_est( seg_cvg )
    [axx, mxi] = max( seg_cvg );
    seq_int8 = int8( mxi-1 );
end

%% f04_fastq_get_read_v07
function [rd_seq, rd_len, Q_av, Qmin, Nmin, n_rds] = f04_fastq_get_read_v09( fp_fa, Min_rd_len, Min_rd_len2, Min_tail_length, n_d_th, minQ, dsf, fp_sr, Nominal_read_length, rdf_en )

n_rds = 0;
min_tail_length = Min_tail_length;
n_line = 0;
% while(1)
    if dsf > 1
        for k = 1:1:dsf-1
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
        end
    end
    cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
    if isempty(cstr{1}) 
        rd_len = -1;
        rd_seq = [];
        Q_av = -1;
        Qmin = -1;
        Nmin = -1;
        % break;
    else
        line_a = cstr{1}{:};
        rd_len = -1;
        rd_seq = [];
        Q_av = -1;
        Qmin = -1;
        Nmin = -1;
        if line_a(1) == '@'

            cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr{1})
                % break;
            else
            line_b = upper( cstr{1}{:} );
            cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr{1})
                % break;
            else
            line_c = cstr{1}{:};
            cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr{1})
                % break;
            else
            line_d = cstr{1}{:};
            n_rds = n_rds + 1;
            if length(line_b) > 2
                if line_b(2) == 'N' || line_b(2) == 'n'
                    sp = 3;
                else
                    sp = 1;
                end
            else
                sp = 1;
            end
            if length(line_b) < 3
                ep = 0;
            else
                if rdf_en == 0
                    ep = length(line_b);
                else
                      ep2 = f04_find_atail_v01( line_d(sp:length(line_b)) );
                      % [sp2, ep2] = f04_find_atail_v02( line_b(sp:length(line_b)) );
                      ep = sp + ep2 - 1; 
                      % sp = sp + sp2 - 1;
                end
            end
            NTstr = ( line_b(sp:ep) );
            tail_len = 0;
            % bin = 0;
            len = length(NTstr);
            if len < Min_rd_len
                bin = 3;
            else
                if Min_tail_length > 0
                    [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( NTstr, Min_tail_length, n_d_th, Nominal_read_length );
                    tail_len = max(tail_len_f, tail_len_b);
                    if tail_len > 0 && tail_len_f < tail_len_b
                        tail_len = -tail_len;
                    end
                    if tail_len == 0
                        bin = 0;
                    else
                        t_tmp = 0;
                        if tail_len_f > 0
                            if sum( NTstr ~= 'A' ) <= Nominal_read_length*0.2 %Min_rd_len2
                                t_tmp = 1;
                            end
                        else
                            if sum( NTstr ~= 'T' ) <= Nominal_read_length*0.2 %Min_rd_len2
                                t_tmp = 1;
                            end
                        end
                        if t_tmp > 0
                            bin = 3;
                        else
                            if len-abs(tail_len) < Min_rd_len2
                                bin = 2;
                            else
                                bin = 1;
                            end
                        end
                    end
                else
                    bin = 0;
                end
            end
            tl = tail_len;
            
            Q_av = -1;
            Qmin = -1;
            Nmin = -1;

            switch( bin )
                case 3,
                    rd_len = 0;
                    rd_seq = [];
                    if fp_sr > 0
                        fprintf(fp_sr, '%s', line_a );
                        fprintf(fp_sr, '%s', line_b );
                        fprintf(fp_sr, '%s', line_c );
                        fprintf(fp_sr, '%s', line_d );
                    end
                case 2,
                    rd_len = 0;
                    rd_seq = [];
                    if fp_sr > 0
                        fprintf(fp_sr, '%s', line_a );
                        fprintf(fp_sr, '%s', line_b );
                        fprintf(fp_sr, '%s', line_c );
                        fprintf(fp_sr, '%s', line_d );
                    end
                case 1,
                    [rd_seq_tmp, rd_len, n_undef] = sub_NTstr2NumSeq( NTstr );  
                    if n_undef > 0
                        rd_len = 0;
                        rd_seq = [];
                    else
                        cut_len = abs(tl) - min_tail_length;
                        if tl > 0
                            rd_len = rd_len - cut_len;
                            rd_seq = rd_seq_tmp(1:rd_len);  
                            qv = uint8( line_d(sp:ep) ) - 33;
                            Q_av = mean( qv );
                            Qmin = min( qv );
                            Nmin = sum( minQ > qv );
                        else
                            rd_seq = 3 - rd_seq_tmp(rd_len:-1:cut_len+1);            
                            rd_len = rd_len - cut_len;
                            qv = uint8( line_d(sp:ep) ) - 33;
                            Q_av = mean( qv );
                            Qmin = min( qv );
                            % Nmin = sum( sign( max(minQ - qv,0) ) );
                            Nmin = sum( minQ > qv );
                         end
                        n_line = n_line+1;
                    end
                otherwise,
                    [rd_seq, rd_len, n_undef] = sub_NTstr2NumSeq( NTstr );  
                    if n_undef > 0
                        rd_len = 0;
                        rd_seq = [];
                    else
                        qv = uint8( line_d(sp:ep) ) - 33;
                        Q_av = mean( qv );
                        Qmin = min( qv );
                        % Nmin = sum( sign( max(minQ - qv,0) ) );
                        Nmin = sum( minQ > qv );
                        n_line = n_line+1;
                    end
            end
            end
            end
            end
        else
            
        end
    end
    if n_line > 0
        % break;
    end
% end
end % sub_fasta_get_read_v05( fp_fa, cfg )

%% f04_fastq_get_read_v07
function [rd_seq, rd_len, Cdep, Q_av, Qmin, Nmin, k_time] = f04_fasta_get_read_v09( fp_fa, Min_rd_len, dsf, flag, min_tail_length, n_d_th, Nominal_read_length, rdf_en )

rd_len = -1;
rd_seq = [];
Cdep = -1;
Q_av = -1;
Qmin = -1;
Nmin = 0;
k_time = 0;

n_line = 0;
% while(1)
    if flag == 0 && dsf > 1
        for k = 1:1:dsf-1
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
        end
    end
    cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
    
    if isempty(cstr{1})
        rd_len = -1;
        rd_seq = [];
        Cdep = -1;
        Q_av = -1;
        Qmin = -1;
        % break;
    else
       line_a = cstr{1}{:};
       if line_a(1) == '>'
            if flag == 0
                Cdep = 0;
                k_time = 0;
                Q_av = 0;
                Qmin = 0;
                Nmin = 0;
            else
                ptr = 2;
                val = sscanf( line_a(ptr:end), '%f' );
                Cdep = val(1);
                k_time = 0;
                Q_av = 0;
                Qmin = 0;
                Nmin = 0;
                if length(val) > 1
                    k_time = val(2);
                    if length(val) > 2
                        Q_av = val(3);
                        if length(val) > 3
                            Qmin = val(4);
                            if length(val) > 4
                                Nmin = val(5);
                            end
                        end
                    end
                end
            end
            
            cstr = textscan(fp_fa,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr{1})
                rd_len = -1;
                rd_seq = [];
                Cdep = -1;
                Q_av = -1;
                Qmin = -1;
                % break;   
            else
                line_b = upper( cstr{1}{:} );
                if flag == 0
                    if length(line_b) > 2
                        if line_b(2) == 'N'
                            sp = 3;
                        else
                            sp = 1;
                        end
                    else
                        sp = 1;
                    end
                    if rdf_en == 0
                        NTstr = line_b(sp:end);
                    else
                        [sp2, ep2] = f04_find_atail_v02( line_b(sp:length(line_b)) );
                        ep = sp + ep2 - 1; 
                        sp = sp + sp2 - 1;
                        NTstr = line_b(sp:ep);
                    end
                else
                    NTstr = line_b;
                end
                if length(NTstr) < Min_rd_len
                    rd_len = 0;
                    rd_seq = [];
                else
                    % tail_len = 0;
                    if flag == 0 && min_tail_length > 0
                        [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( NTstr, min_tail_length, n_d_th, Nominal_read_length );
                        tail_len = max(tail_len_f, tail_len_b);
                        if tail_len > 0 && tail_len_f < tail_len_b
                            tail_len = -tail_len;
                        end
                        if tail_len == 0
                            bin = 0;
                        else
                            t_tmp = 0;
                            if tail_len_f > 0
                                if sum( NTstr ~= 'A' ) <= Nominal_read_length*0.2 %Min_rd_len2
                                    t_tmp = 1;
                                end
                            else
                                if sum( NTstr ~= 'T' ) <= Nominal_read_length*0.2 %Min_rd_len2
                                    t_tmp = 1;
                                end
                            end
                            if t_tmp > 0
                                bin = 3;
                            else
                                if length(NTstr)-abs(tail_len) < 30 %Min_rd_len2
                                    bin = 2;
                                else
                                    bin = 1;
                                end
                            end
                        end
                    else
                        tail_len = 0;
                        bin = 0;
                    end
                    tl = tail_len;
                    switch( bin )
                        case 3,
                            rd_len = 0;
                            rd_seq = [];
                        case 2,
                            rd_len = 0;
                            rd_seq = [];
                        case 1,
                            [rd_seq_tmp, rd_len, n_undef] = sub_NTstr2NumSeq( NTstr );  
                            if n_undef > 0
                                rd_len = 0;
                                rd_seq = [];
                            else
                                cut_len = abs(tl) - min_tail_length;
                                if tl > 0
                                    rd_len = rd_len - cut_len;
                                    rd_seq = rd_seq_tmp(1:rd_len);  
                                else
                                    rd_seq = 3 - rd_seq_tmp(rd_len:-1:cut_len+1);            
                                    rd_len = rd_len - cut_len;
                                 end
                                n_line = n_line+1;
                            end
                        otherwise,
                            [rd_seq, rd_len, n_undef] = sub_NTstr2NumSeq( NTstr );  
                            if n_undef > 0
                                rd_len = 0;
                                rd_seq = [];
                            else
                                n_line = n_line+1;
                            end
                    end
                end
            end
       else
            
        end
    end
    if n_line > 0
        % break;
    end
% end
end % sub_fasta_get_read_v05( fp_fa, cfg )

%% f04_check_polyA_tail_Num_v03

function [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( seq, Min_tail_length, d_threshold, n_rd_len )

    persistent d_th;
    if isempty( d_th )
        d_th = ceil( (1:1:n_rd_len).*d_threshold );
        for k = 1:1:Min_tail_length-1
            d_th(k) = d_th(Min_tail_length);
        end
    end
    Len = length(seq);
    
    if Len > Min_tail_length
        ts = sum( abs(sign( seq(1:Min_tail_length) - 'T' )) );
        if ts > d_th(Min_tail_length);
            tail_len_b = 0;
        else
            b_tmp = 0;
            for k = 1:1:Len-Min_tail_length
                if seq(k+Min_tail_length) ~= 'T'
                    ts = ts + 1;
                    if ts > d_th(k+Min_tail_length) 
                        b_tmp = 1;
                        break;
                    end
                end
            end
            if b_tmp == 0
                tail_len_b = Len;
            else
                if k == 1
                    tail_len_b = 0;
                else
                    if seq(k+Min_tail_length-1) == 'T'
                        tail_len_b = k+Min_tail_length-1;
                    else
                        tail_len_b = k+Min_tail_length-2;
                    end
                end
            end
        end    

        b_tmp = 0;
        ts = sum( abs(sign( seq(Len-Min_tail_length+1:Len) - 'A' )) );
        if ts > d_th(Min_tail_length)
            tail_len_f = 0;
        else
            for k = 1:1:Len-Min_tail_length
                if seq(Len-Min_tail_length+1-k) ~= 'A'
                    ts = ts + 1;
                    if ts > d_th(k+Min_tail_length) 
                        b_tmp = 1;
                        break;
                    end
                end
            end
            if b_tmp == 0
                tail_len_f = Len;
            else
                if k == 1
                    tail_len_f = 0;
                else
                    if seq(Len-Min_tail_length-k+2) == 'A'
                        tail_len_f = k+Min_tail_length-1;
                    else
                        tail_len_f = k+Min_tail_length-2;
                    end
                end
            end
        end
    else
        ts1 = sum( abs(sign( seq(1:Len) - 'T' )) );
        ts2 = sum( abs(sign( seq(1:Len) - 'A' )) );
        if ts1 > d_th(Len) && ts2 > d_th(Len)
            tail_len_f = 0;
            tail_len_b = 0;
        else
            if ts1 <= d_th(Len)
                tail_len_f = Len;
                tail_len_b = Len;
            else
                tail_len_f = Len;
                tail_len_b = Len;
            end
        end
    end
end

%% MISC functions

function ep = f04_find_atail_v01( qseq )
    ep = length(qseq);
    b_tmp = 0;
    for k = 0:1:ep-1
        if uint8(qseq(ep-k)) > 35 
            ep = ep - k;
            b_tmp = 1;
            break;
        end
    end
    if b_tmp == 0
        ep = 1;
    end
end

function [sp, ep] = f04_find_atail_v02( qseq )
    if sum( qseq == 'N' ) == 0
        sp = 1;
        ep = length(qseq);
    else
        ep = length(qseq);
        if qseq(1) == 'N' && qseq(ep) == 'N'
            sp = 1;
            ep = 1;
        else
            if qseq(ep) ~= 'N'
                ep = length(qseq);
                sp = 0;
                for k = 1:1:ep-1
                    if qseq(ep-k) == 'N'
                        sp = ep - (k-1);
                        break;
                    end
                end
                if (ep - sp) < (ep/2) 
                    if qseq(1) ~= 'N'
                        sp = 1;
                        for k = 0:1:ep-1
                            if qseq(sp+k) == 'N' 
                                ep = sp + k;
                                break;
                            end
                        end
                    else
                        
                    end
                end
            else % qseq(1) ~= 'N'
                ep = 0;
                sp = 1;
                for k = 1:1:length(qseq)-1
                    if qseq(sp+k) == 'N'
                        ep = sp + (k-1);
                        break;
                    end
                end
                if (ep - sp) < (ep/2) 
                    if qseq(end) ~= 'N'
                        ep = length(qseq);
                        sp = 0;
                        for k = 1:1:ep-1
                            if qseq(ep-k) == 'N' 
                                sp = ep - (k-1);
                                break;
                            end
                        end
                    else
                    end
                end
            end
        end
    end
end

function tsec = clk2sec( c )
tmon = c(2);
tday = tmon*30 + c(3);
thour = tday*24 + c(4);
tmin = thour*60 + c(5);
tsec = tmin*60 + c(6);
end

function clk = sec2clk( isec )
tsec = mod(isec, 60);
imin = (isec - tsec)/60;
tmin = mod(imin, 60);
ihour = (imin - tmin)/60;
thour = mod(ihour, 24);
iday = (ihour - thour)/24;
tday = mod( iday, 30 );
imon = (iday - tday)/30;
tmon = mod(imon, 12);
iyear = (imon - tmon)/12;
clk = [iyear tmon tday thour tmin tsec];
end

function [n_reads, type, fname, ext, rd_length] = est_n_reads( fname_ext )

rd_length = 0;
type = 0;
b_tmp = 0;
for k = length(fname_ext):-1:1
    if fname_ext(k) == '.'
        b_tmp = 1;
        break;
    end
end
if b_tmp == 0 || k == 1
    fname = fname_ext;
    ext = [];
else
    fname = fname_ext(1:k-1);
    ext = fname_ext(k+1:end);
end

fp = fopen( fname_ext );
if fp < 0
    n_reads = -1;
    fprintf('Cannot open file(s) %s, ', fname_ext );
else
    d = dir( fname_ext );
    n_bytes = d.bytes;
    nrt = 1000;
    n_bytes_read = 0;
    
    d_reads = fgetl(fp);
    n_bytes_read = n_bytes_read + length(d_reads);
    if d_reads(1) == '@'
        type = 2;
    else
        if d_reads(1) == '>'
            type = 1;
        else
            type = 0;
        end
    end
    if type > 0
        for k = 2:1:nrt
            d_reads = fgetl(fp);
            n_bytes_read = n_bytes_read + length(d_reads);
            if mod( k-1, 4 ) == 1
                rd_length = max( rd_length, length(d_reads) );
            end
        end
        fclose(fp);
        if fname_ext(end) == 'q'
            n_bytes_per_read = n_bytes_read/(nrt/4);
        else
            n_bytes_per_read = n_bytes_read/(nrt/2);
        end
        n_reads = round( n_bytes/n_bytes_per_read );
    else
        n_reads = 0;
    end
end
end

function [Num_int8, len, n_undef] = sub_NTstr2NumSeq( NTstr )
persistent ctbl;
if isempty( ctbl )
    ctbl = int8( ones(1,256).*(-1) );
    ctbl( uint8('A') ) = 0;
    ctbl( uint8('C') ) = 1;
    ctbl( uint8('G') ) = 2;
    ctbl( uint8('T') ) = 3;
    ctbl( uint8('a') ) = 0;
    ctbl( uint8('c') ) = 1;
    ctbl( uint8('g') ) = 2;
    ctbl( uint8('t') ) = 3;
end
len = length( NTstr );
Num_int8 = ctbl( uint8(NTstr) ); 
n_undef = sum( max( -Num_int8, 0 ) );
end


function [Nmer_seq, len] = sub_NumSeq2NmerSeq( Num_seq, num_N_mer )

persistent Mask;
if isempty(Mask)
    Mask = uint32(3);
    for k = 1:1:num_N_mer-1
        Mask = bitshift( Mask, 2 ) + 3;
    end
end
if length(Num_seq) < num_N_mer
    len = 1;
    Nmer_seq = uint32( 0 ); 
    for k = 1:1:length(Num_seq)
        Nmer_seq = bitshift( Nmer_seq, 2 ) + uint32( Num_seq(k) );
    end
else
    len = length(Num_seq) - num_N_mer +1;
    Nmer_seq = uint32( zeros(len,1) );
    % First NT -> MSB side of Nmer
    % Last NT -> LSB side of Nmer
    Nmer_seq(1) = uint32( 0 ); 
    for k = 1:1:num_N_mer
        Nmer_seq(1) = bitshift( Nmer_seq(1), 2 ) + uint32( Num_seq(k) );
    end
    for k = 1:1:len-1
        Nmer_seq(k+1) = bitand( Mask, bitshift( Nmer_seq(k), 2 ) + uint32(Num_seq(num_N_mer+k)) );
    end
end
end

function NTstr = sub_NumSeq2NTstr( Num_int8 )

persistent Chs;
if isempty(Chs)
    Chs = 'ACGT';
end
NTstr = Chs( Num_int8 + 1 ); 

end
