
function fname_out = trareco_js( config_file, soption ) 

if exist('soption', 'var') == 0
    soption = 0;
end

if ischar(config_file)
    cfg = trareco_get_config( config_file );
    if isnumeric(cfg) 
        return;
    end
else
    return;
end
out_file_prefix = cfg.output_prefix;

%% load Segments 
dstr = sprintf('Junction search ... started %s', char( datestr(now) ) );
fprintf('%s\n', dstr);

% out_file_prefix = cfg.output_prefix;if ischar( cfg.output_dir )
if ~isempty( cfg.output_dir ) && ischar( cfg.output_dir )
    if ~exist( cfg.output_dir, 'dir' )
        fprintf( '   WARNING: The directory you specified not exists. \n' );
        return;
    end
    out_file_prefix = sprintf('%s/%s', cfg.output_dir, out_file_prefix ); 
end
if soption == 0
    fname_sglb = sprintf('%s.cntg', out_file_prefix );
else
    fname_sglb = sprintf('%s.cntg2', out_file_prefix );
end
str_disp = sprintf('   Loading (%s) ', fname_sglb );
fprintf( '%s', str_disp );
% load( fname_sglb );

%% Set parameters

% Overlap_threshold = cfg.connection_threshold; %min_overlap_depth;
Max_overlap_window = cfg.max_overlap_depth;
Max_num_sel = 100;
Min_sel_th = cfg.cntg_sel_threshold;
Min_tail_length = cfg.connection_threshold_short - 12; %cfg.connection_threshold - 10; %ovlp_th_short - 10; %cfg.min_tail_length;
Junction_overlap_backoff = cfg.junction_backoff;
Max_buff_length = 80000;
Min_seg_length = cfg.nominal_read_length - 1; %75;
Run_thresh = cfg.safe_overlap_threshold;
Normalized_Dist_threshold = cfg.norm_dist_threshold; 

b_inter_disp = 0;
num_N_mer = 8;
if cfg.num_nmer_div_js >= 1 && cfg.num_nmer_div_js <= 8
    seg_sel_dsr = round(cfg.num_nmer_div_js);
else
    seg_sel_dsr = 4;
end
% seg_sel_dsr = 3;

d_scale_ov = 1;
d_scale_cn = 1;
% Blk_len_mul = 3;
Blk_len = round(cfg.nominal_read_length*1.5)*4; %200;
Nx_disp = 100;
Num_proc_sim = cfg.proc_unit_size;
cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth*2) cfg.connection_threshold cfg.bool_ss_ind];

% Ave_num_error_normalized = (cfg.norm_dist_threshold)*(Overlap_threshold);
% Seg_select_threshold = max( (Overlap_threshold - (1.2+Ave_num_error_normalized)*num_N_mer), 10 );
% This threshold is to lower computation time
% (Actually, computation time is not as critical as in contig growing phase)
% For better performance, this threshold must be set as small as possible

if soption == 0
    % fname_out = sprintf('%s_js%s_c%dj%dn%d', out_file_prefix, tver, cfg.connection_threshold_short, cfg.connection_threshold, round(cfg.max_n_cntg_per_group) );
    % fname_out = sprintf('%sn%d', out_file_prefix, round(cfg.max_n_cntg_per_group) );
    fname_out = sprintf('%s', out_file_prefix );
    fname_log = sprintf('%s.log', fname_out );
    fp_log = fopen( fname_log, 'a' );

    if fp_log ~= -1
        fprintf(fp_log, '\n%s\n', dstr);
        fprintf(fp_log, '   # Input: %s.cntg\n', out_file_prefix );
        if cfg.b_split_merged == 0
        %     fprintf(fp_log, '(Possibly) Merged Path split' );
        else
        end
        fprintf(fp_log, '\n%s\n', dstr);
        fprintf(fp_log, '   # Input: %s.cntg\n\n', out_file_prefix );
    end
else
    % fname_out = sprintf('%sj%dn%d', out_file_prefix, cfg.connection_threshold, round(cfg.max_n_cntg_per_group) );
    % fname_out = sprintf('%sn%d', out_file_prefix, round(cfg.max_n_cntg_per_group) );
    fname_out = sprintf('%s', out_file_prefix );
    fname_log = sprintf('%s.log', fname_out );
    fp_log = fopen( fname_log, 'a' );
end

%% Set running variables
if soption == 0
    fname_txt = sprintf('%s.cntg', out_file_prefix );
else
    fname_txt = sprintf('%s.cntg2', out_file_prefix );
end
fp_t = fopen( fname_txt, 'rt' );
n_cnt = 0;
max_clen = 0;
while(1)
    linea = fgets(fp_t);
    lineb = fgets(fp_t);
    fgets(fp_t);
    fgets(fp_t);
    fgets(fp_t);
    linef = fgets(fp_t);
    if linef < 0
        break;
    else
        if linea(1) == '>'
            n_cnt = n_cnt + 1;
            clen = length(lineb)-1;
            max_clen = max( max_clen, clen );
        end
    end
    if mod(n_cnt, 50000) == 0
        fprintf('.' );
    end
end
fclose(fp_t);
% n_cnt

n_seg = n_cnt;
% seg_len_max = max_clen;

sg.id = uint32(0);
sg.tlen = uint32(0);
sg.len = uint32(0);
sg.seq = int8( [] );
sg.cvg_dep = uint32( zeros(4,0) );
sg.ave_cvg_dep = uint32(0);
seg = repmat( sg, n_seg*2, 1 );

fp_t = fopen( fname_txt, 'rt' );
n_cnt = 0;
kmod = round(n_seg/5);

for k = 1:1:n_seg
    linea = fgets(fp_t);
    lineb = fgets(fp_t);
    linec = fgets(fp_t);
    lined = fgets(fp_t);
    linee = fgets(fp_t);
    linef = fgets(fp_t);
    if linef < 0
        break;
    else
        if linea(1) == '>'
            n_cnt = n_cnt + 1;
            
            ptr = 1;
            [axx, bxx, cxx, next_idx] = sscanf( linea(ptr:end), '%s', 1 );
            ptr = ptr + next_idx;
            [val, axx, bxx, cxx] = sscanf( linea(ptr:end), '%d', 4 );
            N_reads_valid = val(3);
            N_reads_total = val(4);

            seg(n_cnt).id = uint32( val(2) );
            seg(n_cnt).len = length(lineb)-1;
            % seg(n_cnt).seq = zeros( 1, seg(n_cnt).len + 1000;
            seg(n_cnt).seq = sub_NTstr2NumSeq( lineb(1:end-1) );
            seg(n_cnt).cvg_dep = zeros( 4, seg(n_cnt).len, 'uint32' );
            %for m = 1:1:4
            consen_tmp = sscanf( linec, '%d', [1 seg(n_cnt).len] );
            seg(n_cnt).cvg_dep(1,:) = consen_tmp;
            consen_tmp = sscanf( lined, '%d', [1 seg(n_cnt).len] );
            seg(n_cnt).cvg_dep(2,:) = consen_tmp;
            consen_tmp = sscanf( linee, '%d', [1 seg(n_cnt).len] );
            seg(n_cnt).cvg_dep(3,:) = consen_tmp;
            consen_tmp = sscanf( linef, '%d', [1 seg(n_cnt).len] );
            seg(n_cnt).cvg_dep(4,:) = consen_tmp;
            %end
            seg(n_cnt).ave_cvg_dep = mean( sum( seg(n_cnt).cvg_dep ) );
        end
    end
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
fclose(fp_t);

if soption == 0
    fname_rcm = sprintf('%s.rcmap', out_file_prefix );
else
    fname_rcm = sprintf('%s.rcmap2', out_file_prefix );
end
fp_r = fopen( fname_rcm, 'rt' );
cntg_cnt = fscanf(fp_r,'%d', 1 );
cntg_valid_ind = fscanf(fp_r,'%d', cntg_cnt );
N_reads_total = fscanf(fp_r,'%d', 1 );
read_cntg_map = zeros( N_reads_total, 3 );
for k = 1:1:N_reads_total
    val = fscanf(fp_r,'%d', 2 );
    read_cntg_map(k,1:2) = val';
end
% smt = sum( read_cntg_map(1:N_reads_total,:) );
% if smt(1) == 0 || smt(2) == 0 
if cfg.read_mode == 0 
    R_mode = 0;
else
    R_mode = 1;
    N_reads_total = N_reads_total*2;
end
fclose(fp_r);
fprintf(' RMode: %d, Nrt: %d', R_mode, N_reads_total );

%% Check integrity
seg_cnt = 0;
for n = 1:n_seg
    [axx,nc] = size( seg(n).cvg_dep );
    if nc ~= seg(n).len || seg(n).len ~= length(seg(n).seq)
        seg_cnt = seg_cnt + 1;
    end
end
if seg_cnt > 0
    str = sprintf('\n ERROR A: Check failed ');
    fprintf('%s', str);
end

% str_disp = sprintf(' %d contigs read, MxLen: %d', n_seg, seg_len_max );
% % fprintf( '%s', str_disp );
% disp( str_disp );
% 

%% Option

if soption == 0

%% Remove PolyA tail and set tail indicator 

% str_disp = sprintf('   Cut PolyA tail ... ' );
% fprintf( '%s', str_disp );

b_seg_valid = ones(1,n_seg);
tail_ind = zeros(1,n_seg);
D_threshold_tail = cfg.norm_dist_threshold; 
mx_seg_len_tmp = max( [seg(1:n_seg).len] );
for n = 1:1:n_seg
    [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( sub_NumSeq2NTstr( seg(n).seq ), Min_tail_length, D_threshold_tail, mx_seg_len_tmp );
    max_tail_len = seg(n).len;
    if tail_len_f == 0 && tail_len_b == 0
        % No action
        tail_ind(n) = 0;
    else
        if abs(tail_len_f) == max_tail_len || abs(tail_len_b) == max_tail_len % - num_N_mer
            if b_inter_disp > 0
                str_disp = sprintf('      ERROR: PolyA tail length %d larger than Max value %d', tail_ind(n), max_tail_len );
                disp( str_disp );
            end
            b_seg_valid(n) = 0;
        else
            if tail_len_b > 0 && tail_len_f > 0
                if tail_len_f > tail_len_b
                    seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                    tail_ind(n) = Min_tail_length;
                    seg(n).cvg_dep = seg(n).cvg_dep(:,1:seg(n).len);
                    seg(n).seq = f03_seq_est( seg(n).cvg_dep );
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep ) );
                else
                    seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                    tail_ind(n) = Min_tail_length;
                    seg(n).cvg_dep = seg(n).cvg_dep(4:-1:1,(tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                    seg(n).seq = 3-seg(n).seq((tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep ) );
                end
            else
                if tail_len_b > 0 || tail_len_f > 0
                    if tail_len_f > tail_len_b
                        seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                        tail_ind(n) = Min_tail_length;
                        seg(n).cvg_dep = seg(n).cvg_dep(:,1:seg(n).len);
                        seg(n).seq = f03_seq_est( seg(n).cvg_dep );
                        seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep ) );
                    else
                        seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                        tail_ind(n) = Min_tail_length;
                        seg(n).cvg_dep = seg(n).cvg_dep(4:-1:1,(tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                        seg(n).seq = 3-seg(n).seq((tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                        seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep ) );
                    end
                end
            end
        end
    end
end
% tail_ind = tail_ind./max_tail_len;

seg_len_ave = 0;
seg_len_max = 0;
seg_len_min = 100000;
seg_cvg_ave = 0;
seg_cvg_max = 0;
seg_cvg_min = 100000;
nseg_len_ave = 0;
nseg_len_min = 0;
nseg_cvg_ave = 0;
nseg_cvg_min = 0;

seg_cnt = 0;
mf = 1.6;
for n = 1:n_seg
    if b_seg_valid(n) == 1 && seg(n).len > num_N_mer
        seg_cnt = seg_cnt + 1;
        % if seg_cnt ~= n
            seg(seg_cnt).id = seg(n).id;
            seg(seg_cnt).len = seg(n).len;
            seg(seg_cnt).tlen = round(seg(seg_cnt).len*mf);
            cvg_tmp = seg(n).cvg_dep;
            seg(seg_cnt).cvg_dep = zeros( 4, seg(seg_cnt).tlen, 'uint32' );
            seg(seg_cnt).cvg_dep(:,1:seg(seg_cnt).len) = cvg_tmp;
            seq_tmp = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(seg_cnt).len) );
            seg(seg_cnt).seq = zeros( 1, seg(seg_cnt).tlen, 'int8' );
            seg(seg_cnt).seq(1:seg(seg_cnt).len) = seq_tmp;
            seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(seg_cnt).len) ) );
            tail_ind(seg_cnt) = tail_ind(n);
        % end
        seg_len_ave = seg_len_ave + seg(seg_cnt).len;
        seg_len_max = max( seg_len_max, seg(seg_cnt).len);
        seg_len_min = min( seg_len_min, seg(seg_cnt).len);
        seg_cvg_ave = seg_cvg_ave + seg(seg_cnt).ave_cvg_dep;
        seg_cvg_max = max( seg_cvg_max, seg(seg_cnt).ave_cvg_dep);
        seg_cvg_min = min( seg_cvg_min, seg(seg_cnt).ave_cvg_dep);
    else
        cntg_valid_ind( seg(n).id ) = -2; % discard
        b_seg_valid(n) = 0;
    end
end
n_seg = seg_cnt;
% fprintf( ' Done \n' );

str_disp = sprintf(', %d contigs read, MxLen: %d ', n_seg, seg_len_max );
fprintf( '%s', str_disp );
% disp( str_disp );

seg_len_ave = seg_len_ave/n_seg;
seg_cvg_ave = seg_cvg_ave/n_seg;
for n = 1:n_seg
    nseg_len_ave = nseg_len_ave + double( seg(n).len <= seg_len_ave );
    nseg_len_min = nseg_len_min + double( seg(n).len <= seg_len_min );
    nseg_cvg_ave = nseg_cvg_ave + double( seg(n).ave_cvg_dep <= seg_cvg_ave );
    nseg_cvg_min = nseg_cvg_min + double( seg(n).ave_cvg_dep <= seg_cvg_min );
end

fprintf(fp_log, '      Num of contigs: %d\n', n_seg );
fprintf(fp_log, '      Maximum contig length: %d\n', seg_len_max );
fprintf(fp_log, '      Average contig length: %f, (%f percent below ave.) \n', seg_len_ave, nseg_len_ave*100/n_seg );
fprintf(fp_log, '      Minimum contig length: %d (%f percent) \n', seg_len_min, nseg_len_min*100/n_seg );
fprintf(fp_log, '      Maximum contig coverage: %f\n', seg_cvg_max );
fprintf(fp_log, '      Average contig coverage: %f (%f percent below ave.) \n', seg_cvg_ave, nseg_cvg_ave*100/n_seg );
fprintf(fp_log, '      Minimum contig coverage: %f (%f percent) \n', seg_cvg_min, nseg_cvg_min*100/n_seg );

%% contig combining

str_hdr = sprintf('   Connecting contigs ....' ); %Checking complete overlap ' );
fprintf( '\n' );
fprintf( '%s', str_hdr );
fprintf( fp_log, '%s', str_hdr );

Nchar = 0;
str_disp = sprintf(' %d', n_seg );
fprintf(repmat('\b', 1, Nchar));
fprintf( '%s', str_disp ); 
fprintf( fp_log, '%s', str_disp );

%% remove Segs having cvg depth less than the minimum value
% % if cfg.min_cvg_depth_js > 9
% seg_cnt = 0;
% cd_threshold = 4; % + cfg.min_overlap_depth/cfg.nominal_read_length+0.2; %max( cfg.min_cvg_depth_js-1, 1 );
% for n = 1:n_seg
%     seg(n).seq = f03_seq_est( seg(n).cvg_dep );
%     seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep ));
% 
%     if seg(n).ave_cvg_dep > cd_threshold % && seg(n).len > cfg.nominal_read_length % (cfg.min_cntg_length_js)
%     % if seg(n).ave_cvg_dep > 1 && seg(n).len > (cfg.min_cntg_length_js)
%         seg_cnt = seg_cnt + 1;
%         seg(seg_cnt).len = seg(n).len;
%         seg(seg_cnt).cvg_dep = seg(n).cvg_dep(:,1:seg(n).len);
%         seg(seg_cnt).seq = f03_seq_est( seg(seg_cnt).cvg_dep );
%         seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep ) );
%         tail_ind(seg_cnt) = tail_ind(n);
%     end
% end
% % str_disp = sprintf(' - %d = %d', n_seg-seg_cnt, seg_cnt );
% str_disp = sprintf(' -> %d', seg_cnt );
% fprintf( '%s', str_disp );
% fprintf( fp_log, '%s', str_disp );
% n_seg = seg_cnt;
% % end

%% Connect Contigs
% This step is needed to prevent complicated junctions, 
% which causes wierd graphs in graph construction step.
    
Nchar2 = Nchar;

Max_seg_len_to_connect = Max_buff_length; 
Ovlp_threshold = max( cfg.connection_threshold_short, Min_tail_length + 10 );
Ave_num_error_normalized = ceil(Ovlp_threshold*Normalized_Dist_threshold);
Seg_select_threshold = round( max( (Ovlp_threshold - (num_N_mer-seg_sel_dsr))/seg_sel_dsr - Ave_num_error_normalized*(num_N_mer/seg_sel_dsr), Min_sel_th ) );

tcfg = cfg;
tcfg.connection_threshold = Ovlp_threshold;
tcfg.norm_dist_threshold = tcfg.norm_dist_threshold*d_scale_cn;

% Num_proc_sim = round(cfg.proc_unit_size/2);
N_Loops = ceil( n_seg/ Num_proc_sim );
% Seg_Nmer_map = zeros( n_seg, 4^num_N_mer, 'int8' );
clear Seg_Nmer_map;
snm_size = min( Num_proc_sim, n_seg );
[sorted_len, sorted_idx] = sort( [seg(1:n_seg).len], 'descend' );
clear sorted_len;
% sorted_idx = 1:n_seg;

Nchar = 0;
max_ovlp = 0;
min_ovlp = 1000;
n_connections = 0;
b_seg_valid = ones(1,n_seg);
n_seg_sel = 0;

%% Update cntg_valid_ind
if R_mode > 0
    
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) < 0
        % cntg_gid(k) = 0;
    else
        if cntg_valid_ind(k) == 0 || cntg_valid_ind(k) == k
            cntg_valid_ind(k) = k;
        else
            cid_p = k;
            gid_tmp = 0;
            while(1)
                cid = cntg_valid_ind(cid_p);
                if cntg_valid_ind(cid) < 0
                    gid_tmp = -1;
                    break;
                else
                    if cntg_valid_ind(cid) == 0 || cntg_valid_ind(cid) == cid
                        gid_tmp = cid;
                        break;
                    else
                        cid_p = cid;
                    end
                end
            end
            cid_p = k;
            while(1)
                cid = cntg_valid_ind(cid_p);
                cntg_valid_ind(cid_p) = gid_tmp;
                if cntg_valid_ind(cid) <= 0 || cntg_valid_ind(cid) == cid
                    cntg_valid_ind(cid) = gid_tmp;
                    break;
                else
                    cid_p = cid;
                end
            end
        end
    end
end

%% set cid_to_seg_map using cntg_valid_ind
cid_to_seg_map = zeros( cntg_cnt, 1, 'uint32' );
for k = 1:1:n_seg
    cid_to_seg_map( seg(k).id ) = k;
end
N_error = 0;
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) > 0 % && cntg_valid_ind(k) ~= k
        if cid_to_seg_map( k ) == 0
            cid_to_seg_map( k ) = cid_to_seg_map( cntg_valid_ind(k) );
            if cntg_valid_ind( cntg_valid_ind(k) ) ~= cntg_valid_ind(k)
                N_error = N_error + 1;
            end
        else
            if cid_to_seg_map( k ) ~= cid_to_seg_map( cntg_valid_ind(k) )
                N_error = N_error + 1;
            end
        end
    end
end
if N_error > 0
    fprintf(' ERROR: set cid_to_seg_map .. %d ', N_error );
end

%% Set list of possible junctions - Use read_cntg_map, cntg_valid_ind, cid_to_seg_map
Max_n_junctions_per_cntg = 200;
jpc_tmp.bs = Max_n_junctions_per_cntg;
jpc_tmp.nj = 0;
jpc_tmp.ji = zeros( Max_n_junctions_per_cntg, 1, 'uint32' );
jpc = repmat( jpc_tmp, n_seg, 1 );
n_pj = 0;

% if R_mode > 0
Nrt = round( N_reads_total/(R_mode+1) );
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 && read_cntg_map(k,2) > 0
        if read_cntg_map(k,1) ~= read_cntg_map(k,2)
            cntg1 = cntg_valid_ind( read_cntg_map(k,1) );
            cntg2 = cntg_valid_ind( read_cntg_map(k,2) );
            if cntg1 == 0
                cntg1 = read_cntg_map(k,1);
            end
            if cntg2 == 0
                cntg2 = read_cntg_map(k,2);
            end
            if cntg1 > 0 && cntg2 > 0
                if cntg1 ~= cntg2
                    for m = 1:1:2
                        if m == 1
                            cid1 = cid_to_seg_map( cntg1 );
                            cid2 = cid_to_seg_map( cntg2 );
                        else
                            cid2 = cid_to_seg_map( cntg1 );
                            cid1 = cid_to_seg_map( cntg2 );
                        end
                        b_tmp = 1;
                        if jpc( cid1 ).nj > 0
                            jmatch = find( jpc( cid1 ).ji(1:jpc( cid1 ).nj) == cid2, 1 );
                            if ~isempty( jmatch )
                                b_tmp = 0;
                            end
                        end
                        if b_tmp > 0
                            jpc( cid1 ).nj = jpc( cid1 ).nj + 1;
                            if jpc( cid1 ).nj > jpc( cid1 ).bs
                                jpc( cid1 ).ji = [jpc( cid1 ).ji; zeros( Max_n_junctions_per_cntg, 1, 'uint32' ) ];
                                jpc( cid1 ).bs = jpc( cid1 ).bs + Max_n_junctions_per_cntg;
                            end
                            jpc( cid1 ).ji( jpc( cid1 ).nj ) = cid2;
                            if m == 1
                                n_pj = n_pj +1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% pre-connection

    cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth*1.5) cfg.connection_threshold cfg.bool_ss_ind];
    s_idxs = find( [jpc(1:end).nj] > 0 );
    nn_cnt = 0;
    for n = 1:1:length(s_idxs) % n_seg
        divisor_idx = s_idxs(n); 
        
        SLen = seg(divisor_idx).len;
        if SLen <= Max_seg_len_to_connect && b_seg_valid(divisor_idx) > 0
            
            sel_seg_idx = jpc(divisor_idx).ji(1:jpc(divisor_idx).nj);
            
            % [sorted_match, sorted_idx2] = sort( Nmer_Matches(sel_seg_idx), 'descend' );
            depth_all  = zeros( length(sel_seg_idx)*4, 1 );
            c_mode_all  = zeros( length(sel_seg_idx)*4, 1 );
            dist_all  = ones( length(sel_seg_idx)*4, 1 ).*10000;
            n_seg_sel = 0.99*n_seg_sel + 0.01*length(sel_seg_idx);
                        
            for k = 1:1:length(sel_seg_idx)*4 %n_seg 
                idx_t = ceil( k/4 ) ;
                dividend_idx = sel_seg_idx( idx_t );
                t_mode = mod(k-1,4)+1;
                if (divisor_idx ~= dividend_idx) && (b_seg_valid(dividend_idx)) > 0 % && (seg(dividend_idx).len+SLen <= Max_seg_len_to_connect-cfg.connection_threshold_short) 
                    % [depth_all(k), c_mode_all(k), dist] = sub_tcn( seg(dividend_idx).seq, seg(divisor_idx).seq, cfg_val, t_mode );
                    % cfg_val = [cfg.norm_dist_threshold 16 round(cfg.max_overlap_depth*3) cfg.connection_threshold cfg.bool_ss_ind];
                    [depth_all(k), c_mode_all(k), dist_all(k)] = sub_tcn( ...
                        seg(dividend_idx).seq(1:seg(dividend_idx).len), seg(divisor_idx).seq(1:seg(divisor_idx).len), cfg_val, t_mode );
                    
                    if dist_all(k) == 0 % depth_all(k) > Run_thresh
                        break;
                    end
                end
            end
            [depth, midx] = min( dist_all.*cfg.nominal_read_length - depth_all ); %max( depth_all - dist_all );
            depth = depth_all(midx);
            c_mode = c_mode_all(midx);
            if c_mode ~= 0
                % dividend_idx = sorted_idx( idxs( sel_seg_idx2(sel_seg_idx(midx)) )+sp-1 );
                idx_t = ceil( midx/4 ) ;
                dividend_idx = sel_seg_idx( idx_t );
                n_connections = n_connections + 1;

                if max_ovlp < depth
                    max_ovlp = depth;
                end
                if min_ovlp > depth
                    min_ovlp = depth;
                end
                nn_cnt = nn_cnt +1;
                o_dep = depth;
                idx_r = dividend_idx;
                idx_o = divisor_idx;
                if c_mode < 0
                    seg_idx = idx_o;
                    seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) = seg(seg_idx).cvg_dep(4:-1:1,seg(seg_idx).len:-1:1);
                    seg(seg_idx).seq(1:seg(seg_idx).len) = f03_seq_est( seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) );
                    seg(seg_idx).ave_cvg_dep = mean( sum( seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) ) );
                end
                if abs(c_mode) == 1
                    seg(idx_r).cvg_dep(:,1:o_dep) = ...
                       seg(idx_r).cvg_dep(:,1:o_dep) + seg(idx_o).cvg_dep(:,double(seg(idx_o).len)-o_dep+1:double(seg(idx_o).len));
                    seg_len_new = double(seg(idx_r).len + seg(idx_o).len) - o_dep;
                    if seg_len_new > seg(idx_r).tlen
                        seg(idx_r).tlen = round(seg(idx_r).tlen*mf);
                        cvg_tmp = seg(idx_r).cvg_dep(:,1:seg(idx_r).len);
                        seg(idx_r).cvg_dep = zeros( 4, seg(idx_r).tlen, 'uint32' );
                        seg(idx_r).cvg_dep(:,1:seg(idx_r).len) = cvg_tmp;
                        seq_tmp = f03_seq_est( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) );
                        seg(idx_r).seq = zeros( 1, seg(idx_r).tlen, 'int8' );
                        seg(idx_r).seq(1:seg(idx_r).len) = seq_tmp;
                        seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) ) );
                        tail_ind(idx_r) = tail_ind(idx_r);
                    end
                    seg(idx_r).cvg_dep(:,1:seg_len_new) = [seg(idx_o).cvg_dep(:,1:double(seg(idx_o).len)-o_dep) seg(idx_r).cvg_dep(:,1:seg(idx_r).len)];
                    seg(idx_r).seq(1:seg_len_new) = [seg(idx_o).seq(1:double(seg(idx_o).len)-o_dep) seg(idx_r).seq(1:seg(idx_r).len)];
                    seg(idx_r).len = seg_len_new;
                    seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg_len_new) ) );
                else
                    seg(idx_r).cvg_dep(:,double(seg(idx_r).len)-o_dep+1:double(seg(idx_r).len)) = ...
                       seg(idx_r).cvg_dep(:,double(seg(idx_r).len)-o_dep+1:double(seg(idx_r).len)) + seg(idx_o).cvg_dep(:,1:o_dep);
                    seg_len_new = double(seg(idx_r).len + seg(idx_o).len) - o_dep;
                    if seg_len_new > seg(idx_r).tlen
                        seg(idx_r).tlen = round(seg(idx_r).tlen*mf);
                        cvg_tmp = seg(idx_r).cvg_dep(:,1:seg(idx_r).len);
                        seg(idx_r).cvg_dep = zeros( 4, seg(idx_r).tlen, 'uint32' );
                        seg(idx_r).cvg_dep(:,1:seg(idx_r).len) = cvg_tmp;
                        seq_tmp = f03_seq_est( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) );
                        seg(idx_r).seq = zeros( 1, seg(idx_r).tlen, 'int8' );
                        seg(idx_r).seq(1:seg(idx_r).len) = seq_tmp;
                        seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) ) );
                        tail_ind(idx_r) = tail_ind(idx_r);
                    end
                    seg(idx_r).cvg_dep(:,1:seg_len_new) = [seg(idx_r).cvg_dep(:,1:seg(idx_r).len) seg(idx_o).cvg_dep(:,o_dep+1:double(seg(idx_o).len))];
                    seg(idx_r).seq(1:seg_len_new) = [seg(idx_r).seq(1:seg(idx_r).len) seg(idx_o).seq(o_dep+1:double(seg(idx_o).len))];
                    seg(idx_r).len = seg_len_new;
                    seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg_len_new) ) );
                end
                                
                b_seg_valid(idx_o) = 0;
                if tail_ind(idx_o) > 0
                    tail_ind(idx_r) = tail_ind(idx_o);
                end
                cntg_valid_ind( seg(idx_o).id ) = seg(idx_r).id;
            end
        end
        % n_step = 5; %round(n_seg/10);
        % if mod(n, Nx_disp ) == 0
            % fprintf( '>' );
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            % Nchar = fprintf( ' %d/%d', n, n_seg ); 
            % Nchar = fprintf( ' %d(%d/%d)', n-(sp+1), lp, N_Loops ); 
            Nchar = fprintf( ' .. %d/%d/%d ', nn_cnt, n, length(s_idxs) ); 
        % end
    end

else
    
%% Connecting Contig start
Nchar = 0;
% if Nchar > 0
cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth*1.5) cfg.connection_threshold cfg.bool_ss_ind];
for lp = N_Loops:-1:1 %1:1:N_Loops
    
    Nchar = Nchar + fprintf( '  Wait ' ); 
    sp = (lp-1)*Num_proc_sim+1;
    if lp < N_Loops
        ep = lp*Num_proc_sim;
    else
        ep = n_seg;
    end
    num_proc_sim = ep-sp+1;

    Seg_Nmer_map = zeros( snm_size*4, 4^num_N_mer, 'int8' );
    for n = 1:num_proc_sim
        idx = sorted_idx(n+sp-1);
        SLen = seg(idx).len;
        OLen = min( Max_overlap_window, seg(idx).len );
        [Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:OLen), num_N_mer );
        % Seg_Nmer_map( 4*(n-1)+2,: ) = 0;
        Seg_Nmer_map( 4*(n-1)+2, Nmer1(1:sl1)+1 ) = 1;
        [Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(OLen:-1:1), num_N_mer );
        % Seg_Nmer_map( 4*(n-1)+3,: ) = 0;
        Seg_Nmer_map( 4*(n-1)+3, Nmer2(1:sl2)+1 ) = 1;
        [Nmer3, sl3] = sub_NumSeq2NmerSeq( seg(idx).seq(SLen-OLen+1:SLen), num_N_mer );
        % Seg_Nmer_map( 4*(n-1)+1,: ) = 0;
        Seg_Nmer_map( 4*(n-1)+1, Nmer3(1:sl3)+1 ) = 1;
        [Nmer4, sl4] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(SLen:-1:SLen-OLen+1), num_N_mer );
        % Seg_Nmer_map( 4*(n-1)+4,: ) = 0;
        Seg_Nmer_map( 4*(n-1)+4, Nmer4(1:sl4)+1 ) = 1;
        if mod(n,2000) == 0
            Nchar = Nchar + fprintf( '.' );
        end
    end

    % Connection discovery
    for n = n_seg:-1:sp+1
    % for n = n_seg:-1:2 %1:1:n_seg
        divisor_idx = sorted_idx(n); 
        SLen = seg(divisor_idx).len;
        if SLen <= Max_seg_len_to_connect && b_seg_valid(divisor_idx) > 0
            
            OLen = min( Max_overlap_window, seg(divisor_idx).len );
            [Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(divisor_idx).seq(1:OLen), num_N_mer );
            [Nmer3, sl3] = sub_NumSeq2NmerSeq( seg(divisor_idx).seq(SLen-OLen+1:SLen), num_N_mer );

            Match_threshold = Seg_select_threshold; %round(Seg_select_threshold/seg_sel_dsr);
            n_cmp = min( ep, n-1 ) - sp+1;
            idxs1 = 1:2:n_cmp*4; %find( b_seg_valid( sorted_idx(1:n_cmp) ) > 0 );
            idxs2 = 2:2:n_cmp*4; %find( b_seg_valid( sorted_idx(1:n_cmp) ) > 0 );
            Nmer_Matches1 = sum( Seg_Nmer_map( idxs1, Nmer1(1:seg_sel_dsr:sl1)+1 ), 2 );
            Nmer_Matches2 = sum( Seg_Nmer_map( idxs2, Nmer3(1:seg_sel_dsr:sl3)+1 ), 2 );
            idxs = 1:1:n_cmp*4;
            Nmer_Matches = zeros( 1, n_cmp*4 );
            Nmer_Matches( idxs1 ) = Nmer_Matches1;
            Nmer_Matches( idxs2 ) = Nmer_Matches2;
            
            sel_seg_idx = find( Nmer_Matches >= Match_threshold );
            % [sorted_match, sorted_idx2] = sort( Nmer_Matches(sel_seg_idx), 'descend' );
            depth_all  = zeros( length(sel_seg_idx), 1 );
            c_mode_all  = zeros( length(sel_seg_idx), 1 );
            dist_all  = ones( length(sel_seg_idx), 1 ).*10000;
            n_seg_sel = 0.99*n_seg_sel + 0.01*length(sel_seg_idx);
                        
            for k = 1:1:length(sel_seg_idx) %n_seg 
                idx_t = floor( (idxs(sel_seg_idx(k))-1)/4 );
                dividend_idx = sorted_idx( idx_t + sp );
                t_mode = idxs( sel_seg_idx(k) ) - 4*idx_t;
                if (divisor_idx ~= dividend_idx) && (b_seg_valid(dividend_idx)) > 0 && (seg(dividend_idx).len+SLen <= Max_seg_len_to_connect-cfg.connection_threshold_short) 
                    % [depth_all(k), c_mode_all(k), dist] = sub_tcn( seg(dividend_idx).seq, seg(divisor_idx).seq, cfg_val, t_mode );
                    [depth_all(k), c_mode_all(k), dist_all(k)] = sub_tcn( ...
                        seg(dividend_idx).seq(1:seg(dividend_idx).len), seg(divisor_idx).seq(1:seg(divisor_idx).len), cfg_val, t_mode );
                    
                    if dist_all(k) == 0 % depth_all(k) > Run_thresh
                        break;
                    end
                end
            end
            [depth, midx] = min( dist_all.*cfg.nominal_read_length - depth_all ); %max( depth_all - dist_all );
            depth = depth_all(midx);
            c_mode = c_mode_all(midx);
            if c_mode ~= 0
                % dividend_idx = sorted_idx( idxs( sel_seg_idx2(sel_seg_idx(midx)) )+sp-1 );
                idx_t = floor( (idxs(sel_seg_idx(midx))-1)/4 );
                dividend_idx = sorted_idx( idx_t + sp );
                n_connections = n_connections + 1;

                if max_ovlp < depth
                    max_ovlp = depth;
                end
                if min_ovlp > depth
                    min_ovlp = depth;
                end

                o_dep = depth;
                idx_r = dividend_idx;
                idx_o = divisor_idx;
                if c_mode < 0
                    seg_idx = idx_o;
                    seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) = seg(seg_idx).cvg_dep(4:-1:1,seg(seg_idx).len:-1:1);
                    seg(seg_idx).seq(1:seg(seg_idx).len) = f03_seq_est( seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) );
                    seg(seg_idx).ave_cvg_dep = mean( sum( seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) ) );
                end
                if abs(c_mode) == 1
                    seg(idx_r).cvg_dep(:,1:o_dep) = ...
                       seg(idx_r).cvg_dep(:,1:o_dep) + seg(idx_o).cvg_dep(:,double(seg(idx_o).len)-o_dep+1:double(seg(idx_o).len));
                    seg_len_new = double(seg(idx_r).len + seg(idx_o).len) - o_dep;
                    if seg_len_new > seg(idx_r).tlen
                        seg(idx_r).tlen = round(seg(idx_r).tlen*mf);
                        cvg_tmp = seg(idx_r).cvg_dep(:,1:seg(idx_r).len);
                        seg(idx_r).cvg_dep = zeros( 4, seg(idx_r).tlen, 'uint32' );
                        seg(idx_r).cvg_dep(:,1:seg(idx_r).len) = cvg_tmp;
                        seq_tmp = f03_seq_est( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) );
                        seg(idx_r).seq = zeros( 1, seg(idx_r).tlen, 'int8' );
                        seg(idx_r).seq(1:seg(idx_r).len) = seq_tmp;
                        seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) ) );
                        tail_ind(idx_r) = tail_ind(idx_r);
                    end
                    seg(idx_r).cvg_dep(:,1:seg_len_new) = [seg(idx_o).cvg_dep(:,1:double(seg(idx_o).len)-o_dep) seg(idx_r).cvg_dep(:,1:seg(idx_r).len)];
                    seg(idx_r).seq(1:seg_len_new) = [seg(idx_o).seq(1:double(seg(idx_o).len)-o_dep) seg(idx_r).seq(1:seg(idx_r).len)];
                    seg(idx_r).len = seg_len_new;
                    seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg_len_new) ) );
                else
                    seg(idx_r).cvg_dep(:,double(seg(idx_r).len)-o_dep+1:double(seg(idx_r).len)) = ...
                       seg(idx_r).cvg_dep(:,double(seg(idx_r).len)-o_dep+1:double(seg(idx_r).len)) + seg(idx_o).cvg_dep(:,1:o_dep);
                    seg_len_new = double(seg(idx_r).len + seg(idx_o).len) - o_dep;
                    if seg_len_new > seg(idx_r).tlen
                        seg(idx_r).tlen = round(seg(idx_r).tlen*mf);
                        cvg_tmp = seg(idx_r).cvg_dep(:,1:seg(idx_r).len);
                        seg(idx_r).cvg_dep = zeros( 4, seg(idx_r).tlen, 'uint32' );
                        seg(idx_r).cvg_dep(:,1:seg(idx_r).len) = cvg_tmp;
                        seq_tmp = f03_seq_est( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) );
                        seg(idx_r).seq = zeros( 1, seg(idx_r).tlen, 'int8' );
                        seg(idx_r).seq(1:seg(idx_r).len) = seq_tmp;
                        seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len) ) );
                        tail_ind(idx_r) = tail_ind(idx_r);
                    end
                    seg(idx_r).cvg_dep(:,1:seg_len_new) = [seg(idx_r).cvg_dep(:,1:seg(idx_r).len) seg(idx_o).cvg_dep(:,o_dep+1:double(seg(idx_o).len))];
                    seg(idx_r).seq(1:seg_len_new) = [seg(idx_r).seq(1:seg(idx_r).len) seg(idx_o).seq(o_dep+1:double(seg(idx_o).len))];
                    seg(idx_r).len = seg_len_new;
                    seg(idx_r).ave_cvg_dep = mean( sum( seg(idx_r).cvg_dep(:,1:seg_len_new) ) );
                end
                slen = seg(idx_r).len;
                olen = min( Max_overlap_window, seg(idx_r).len );
                
                idx_n = idx_t*4;
                [Nmer_t1, sl1] = sub_NumSeq2NmerSeq( seg(idx_r).seq(1:olen), num_N_mer );
                Seg_Nmer_map( idx_n+2, : ) = 0;
                Seg_Nmer_map( idx_n+2, Nmer_t1(1:sl1)+1 ) = 1;
                [Nmer_t2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx_r).seq(olen:-1:1), num_N_mer );
                Seg_Nmer_map( idx_n+3, : ) = 0;
                Seg_Nmer_map( idx_n+3, Nmer_t2(1:sl2)+1 ) = 1;
                [Nmer_t3, sl3] = sub_NumSeq2NmerSeq( seg(idx_r).seq(slen-olen+1:slen), num_N_mer );
                Seg_Nmer_map( idx_n+1, : ) = 0;
                Seg_Nmer_map( idx_n+1, Nmer_t3(1:sl3)+1 ) = 1;
                [Nmer_t4, sl4] = sub_NumSeq2NmerSeq( 3-seg(idx_r).seq(slen:-1:slen-olen+1), num_N_mer );
                Seg_Nmer_map( idx_n+4, : ) = 0;
                Seg_Nmer_map( idx_n+4, Nmer_t4(1:sl4)+1 ) = 1;
                
                b_seg_valid(idx_o) = 0;
                if tail_ind(idx_o) > 0
                    tail_ind(idx_r) = tail_ind(idx_o);
                end
                cntg_valid_ind( seg(idx_o).id ) = seg(idx_r).id;
            end
        end
        % n_step = 5; %round(n_seg/10);
        if mod(n-(sp+1), Nx_disp ) == 0
            % fprintf( '>' );
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            % Nchar = fprintf( ' %d/%d', n, n_seg ); 
            % Nchar = fprintf( ' %d(%d/%d)', n-(sp+1), lp, N_Loops ); 
            Nchar = fprintf( ' .. %d/%d(%d/%d)', n-(sp+1), round(n_seg_sel), N_Loops-lp+1, N_Loops ); 
        end
    end
end 
end

if sum( b_seg_valid(1:n_seg) ) ~= n_seg %&& n_depleted == 0
    seg_cnt = 0;
    for n = 1:n_seg
        if b_seg_valid(n) > 0
            seg_cnt = seg_cnt + 1;
            seg(seg_cnt).id = seg(n).id;
            seg(seg_cnt).len = seg(n).len;
            seg(seg_cnt).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
            seg(seg_cnt).seq(1:seg(n).len) = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(n).len) );
            seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(n).len) ) );
            tail_ind(seg_cnt) = tail_ind(n);
        end
    end
    n_seg = seg_cnt;
    fprintf(repmat('\b', 1, Nchar+Nchar2));
    str_disp = sprintf(' -> %d (ML %d/%d, MO %d/%d)', n_seg, max([seg(1:n_seg).len]), min([seg(1:n_seg).len]), min_ovlp, max_ovlp );
    fprintf( '%s', str_disp );
    fprintf( fp_log, '%s', str_disp );
else
    fprintf(repmat('\b', 1, Nchar+Nchar2));
    str_disp = sprintf(' -> %d (ML %d/%d, MO %d/%d)', n_seg, max([seg(1:n_seg).len]), min([seg(1:n_seg).len]), min_ovlp, max_ovlp );
    fprintf( '%s', str_disp );
    fprintf( fp_log, '%s', str_disp );
    % break;
end

%% Remove PolyA tail and set tail indicator 

% str_disp = sprintf('   Cut PolyA tail ... ' );
% fprintf( '\n%s', str_disp );
Min_tail_length = 9; 

if cfg.b_tail_suppress > 0
    b_seg_valid = ones(1,n_seg);
    tail_ind = zeros(1,n_seg);
    D_threshold_tail = cfg.norm_dist_threshold; 
    mx_seg_len_tmp = max( [seg(1:n_seg).len] );
    for n = 1:1:n_seg
        % tail_len = f04_check_polyA_tail_Num_v03( seg_seq_lib(n,1:seg_len(n)), Min_tail_length );
        [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( sub_NumSeq2NTstr( seg(n).seq(1:seg(n).len) ), Min_tail_length, D_threshold_tail, mx_seg_len_tmp );
        max_tail_len = seg(n).len;
        if tail_len_f == 0 && tail_len_b == 0
            % No action
            tail_ind(n) = 0;
        else
            if abs(tail_len_f) == max_tail_len || abs(tail_len_b) == max_tail_len % - num_N_mer
                if b_inter_disp > 0
                    str_disp = sprintf('      ERROR: PolyA tail length %d larger than Max value %d', tail_ind(n), max_tail_len );
                    disp( str_disp );
                end
                b_seg_valid(n) = 0;
            else
                if tail_len_b > 0 && tail_len_f > 0
                    seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                    tail_ind(n) = Min_tail_length;
                    seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,(tail_len_b - Min_tail_length)+1:(tail_len_b - Min_tail_length)+seg(n).len);
                    seg(n).seq(1:seg(n).len) = seg(n).seq((tail_len_b - Min_tail_length)+1:(tail_len_b - Min_tail_length)+seg(n).len);
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                    
                    seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                    seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
                    seg(n).seq(1:seg(n).len) = seg(n).seq(1:seg(n).len);
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                    tail_ind(n) = Min_tail_length;
                else
                    if tail_len_b > 0 || tail_len_f > 0
                        if tail_len_f > tail_len_b
                            seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                            tail_ind(n) = Min_tail_length;
                            seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
                            seg(n).seq(1:seg(n).len) = f03_seq_est( seg(n).cvg_dep(:,1:seg(n).len) );
                            seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                        else
                            seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                            tail_ind(n) = Min_tail_length;
                            seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(4:-1:1,(tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                            seg(n).seq(1:seg(n).len) = 3-seg(n).seq((tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                            seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                        end
                    end
                end
            end
        end
    end
    % tail_ind = tail_ind./max_tail_len;

    seg_cnt = 0;
    for n = 1:n_seg
        if b_seg_valid(n) == 1 && seg(n).len > num_N_mer
            seg_cnt = seg_cnt + 1;
            %if seg_cnt ~= n
                seg(seg_cnt).len = seg(n).len;
                seg(seg_cnt).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
                seg(seg_cnt).seq(1:seg(n).len) = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(n).len) );
                [axx, nc] = size(seg(n).cvg_dep(:,1:seg(n).len));
                if nc ~= seg(n).len || seg(n).len ~= length(seg(seg_cnt).seq(1:seg(n).len))
                    str = sprintf('ERROR 3: %d ~= %d ~= %d', nc, seg(n).len, length(seg(seg_cnt).seq(1:seg(n).len)) );
                    disp( str );
                end
                seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(n).len) ) );
                tail_ind(seg_cnt) = tail_ind(n);
            %end
        else
            cntg_valid_ind( seg(n).id ) = -2;
        end
    end
    n_seg = seg_cnt;
    % fprintf( ' RPT Done ' );
end

%% Merge contigs
Nchar = 0;
Merge_num = 0;
Merge_dep_max = 0;
Merge_dep_ave = 0;

num_N_mer = 8;
seg_sel_dsr = 4;

[sorted_len, seg_idx_sorted] = sort( [seg(1:n_seg).len], 'descend' );
clear sorted_len;
Num_proc_sim = cfg.proc_unit_size*4;
N_Loops = ceil( n_seg/ Num_proc_sim );
snm_size = min( Num_proc_sim, n_seg );

n_seg_merged = 0;
b_seg_valid = ones(1,n_seg);
n_step = Nx_disp;
for lp = 1:1:N_Loops
    
    sp = (lp-1)*Num_proc_sim+1;
    if lp < N_Loops
        ep = lp*Num_proc_sim;
    else
        ep = n_seg;
    end
    num_proc_sim = ep-sp+1;

    Seg_Nmer_map = zeros( snm_size, 4^num_N_mer, 'int8' );
    for n = 1:num_proc_sim %n_seg
        idx = seg_idx_sorted(sp+n-1);
        [seg_seq_Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:seg(idx).len), num_N_mer );
        [seg_seq_Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(seg(idx).len:-1:1), num_N_mer );
        Seg_Nmer_map( n, [seg_seq_Nmer1(1:sl1) seg_seq_Nmer2(1:sl2)]+1 ) = 1;
        % if mod(n, n_step) == 0
        %     fprintf('.' );
        %     % Nchar = Nchar + 1; % fprintf( '%d/%d', (n_seg - n), n_seg ); 
        % end
    end
    for n = n_seg:-1:sp+1
        divisor_idx = seg_idx_sorted(n);
        if b_seg_valid(divisor_idx) > 0
            SLen = seg(divisor_idx).len;
            [Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(divisor_idx).seq(1:SLen), num_N_mer );

            Ave_num_error_normalized = ceil( seg(divisor_idx).len*Normalized_Dist_threshold );
            Match_threshold = round( max( (seg(divisor_idx).len - (num_N_mer-seg_sel_dsr))/seg_sel_dsr - Ave_num_error_normalized*(num_N_mer/seg_sel_dsr), Min_sel_th ) );
            n_cmp = min( ep, n-1 ) - sp+1;
            Nmer_Matches = sum( Seg_Nmer_map( 1:n_cmp, Nmer1(1:seg_sel_dsr:sl1)+1 ), 2 );
            sel_seg_idx2 = find( Nmer_Matches >= Match_threshold );
            sel_seg_idx = sel_seg_idx2( b_seg_valid(seg_idx_sorted(sel_seg_idx2+sp-1)) > 0 );
            if seg(divisor_idx).len < cfg.short_seg_threshold
                Normalized_D_threshold = d_scale_ov*cfg.norm_dist_threshold;
            else
                Normalized_D_threshold = cfg.norm_dist_threshold;
            end
            for k = 1:1:length(sel_seg_idx)
                dividend_idx = seg_idx_sorted(sel_seg_idx(k)+sp-1);
                [pos, ref_seq_idx, dst] = sub_tco( seg(dividend_idx).seq(1:seg(dividend_idx).len), ...
                    seg(divisor_idx).seq(1:seg(divisor_idx).len), Normalized_D_threshold, cfg.bool_ss_ind );
                if pos ~= 0
                    if ref_seq_idx == 0
                        merge_dep = seg(divisor_idx).ave_cvg_dep;
                        Merge_dep_ave = Merge_dep_ave + merge_dep;
                        Merge_num = Merge_num + 1;
                        if Merge_dep_max < merge_dep
                            Merge_dep_max = merge_dep;
                        end
                        if tail_ind(dividend_idx) > 0
                            if tail_ind(divisor_idx) > 0
                                if pos > 0 
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(:,1:seg(divisor_idx).len);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Tail/Len %d/%d merged to Tail/Len %d/%d - %d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len, dst );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                else
                                    pos = -pos;
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(4:-1:1,seg(divisor_idx).len:-1:1);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('complement of Seg/Len %d/%d merged to Tail/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                end
                            else % tail_ind(divisor_idx) == 0
                                if pos > 0 
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(:,1:seg(divisor_idx).len);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Seg/Len %d/%d merged to Tail/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                else
                                    pos = -pos;
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(4:-1:1,seg(divisor_idx).len:-1:1);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('complement of Seg/Len %d/%d merged to Tail/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                               end
                            end
                        else % tail_ind(dividend_idx) == 0
                            if tail_ind(divisor_idx) > 0
                                if pos > 0 
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(:,1:seg(divisor_idx).len);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Tail/Len %d/%d merged to Seg/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                else
                                    pos = -pos;
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(4:-1:1,seg(divisor_idx).len:-1:1);
                                    seg(dividend_idx).cvg_dep(1:4,1:seg(dividend_idx).len) = seg(dividend_idx).cvg_dep(4:-1:1,seg(dividend_idx).len:-1:1);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Tail/Len %d/%d merged to Seg/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                end
                            else % tail_ind(divisor_idx) == 0
                                if pos > 0 
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(:,1:seg(divisor_idx).len);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Seg/Len %d/%d merged to Seg/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                else
                                    pos = -pos;
                                    seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) = ...
                                       seg(dividend_idx).cvg_dep(:,pos:seg(divisor_idx).len+pos-1) + seg(divisor_idx).cvg_dep(4:-1:1,seg(divisor_idx).len:-1:1);
                                    seg(dividend_idx).seq(1:seg(dividend_idx).len) = f03_seq_est( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) );
                                    seg(dividend_idx).ave_cvg_dep = mean( sum( seg(dividend_idx).cvg_dep(:,1:seg(dividend_idx).len) ) );
                                    if b_inter_disp > 0
                                        str_disp = sprintf('Complement of Seg/Len %d/%d merged to Seg/Len %d/%d', ...
                                            divisor_idx, seg(divisor_idx).len, dividend_idx, seg(dividend_idx).len );
                                        fprintf( '\n   %s, ', str_disp );
                                    end
                                end
                            end
                        end
                        b_seg_valid(divisor_idx) = 0;
                        n_seg_merged = n_seg_merged + 1;
                        cntg_valid_ind( seg(divisor_idx).id ) = seg(dividend_idx).id;
                        break;
                    else
                        if b_inter_disp > 0
                            disp( '\n      ERROR: This case cannot be occurred!' );
                        end
                    end
                end
            end
        end
        if mod(n-(sp+1), n_step) == 0
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf( ' .. %d(%d/%d)', n-(sp+1), lp, N_Loops ); 
        end
    end
end

seg_cnt = 0;
for n = 1:n_seg
    if b_seg_valid(n) == 1
        seg_cnt = seg_cnt + 1;
        seg(seg_cnt).id = seg(n).id;
        seg(seg_cnt).len = seg(n).len;
        seg(seg_cnt).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
        [nr, nc] = size(seg(n).cvg_dep(:,1:seg(n).len));
        seg(seg_cnt).seq(1:seg(n).len) = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(n).len) );
        if nc ~= seg(n).len || seg(n).len ~= length(seg(seg_cnt).seq(1:seg(n).len))
            str = sprintf('ERROR 1: %d ~= %d ~= %d, %d', nc, seg(n).len, length(seg(seg_cnt).seq(1:seg(n).len)), nr );
            disp( str );
        end
        seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(n).len) ) );
        tail_ind(seg_cnt) = tail_ind(n);
    end
end
n_seg = seg_cnt; 
str_disp = sprintf(' -> %d', n_seg );
fprintf(repmat('\b', 1, Nchar));
fprintf( '%s', str_disp ); 
fprintf( fp_log, '%s', str_disp );

%% Remove PolyA tail and set tail indicator 

Min_tail_length = 8;  
if cfg.b_tail_suppress > 0
    b_seg_valid = ones(1,n_seg);
    tail_ind = zeros(1,n_seg);
    D_threshold_tail = cfg.norm_dist_threshold; 
    mx_seg_len_tmp = max( [seg(1:n_seg).len] );
    for n = 1:1:n_seg
        [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( sub_NumSeq2NTstr( seg(n).seq(1:seg(n).len) ), Min_tail_length, D_threshold_tail, mx_seg_len_tmp );
        max_tail_len = seg(n).len;
        if tail_len_f == 0 && tail_len_b == 0
            % No action
            tail_ind(n) = 0;
        else
            if abs(tail_len_f) == max_tail_len || abs(tail_len_b) == max_tail_len % - num_N_mer
                if b_inter_disp > 0
                    str_disp = sprintf('      ERROR: PolyA tail length %d larger than Max value %d', tail_ind(n), max_tail_len );
                    disp( str_disp );
                end
                b_seg_valid(n) = 0;
            else
                if tail_len_b > 0 && tail_len_f > 0
                    seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                    tail_ind(n) = Min_tail_length;
                    seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,(tail_len_b - Min_tail_length)+1:(tail_len_b - Min_tail_length)+seg(n).len);
                    seg(n).seq(1:seg(n).len) = seg(n).seq((tail_len_b - Min_tail_length)+1:(tail_len_b - Min_tail_length)+seg(n).len);
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                    
                    seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                    seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
                    seg(n).seq(1:seg(n).len) = seg(n).seq(1:seg(n).len);
                    seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                    tail_ind(n) = Min_tail_length;
                else
                    if tail_len_b > 0 || tail_len_f > 0
                        if tail_len_f > tail_len_b
                            seg(n).len = seg(n).len - (tail_len_f - Min_tail_length);
                            tail_ind(n) = Min_tail_length;
                            seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
                            seg(n).seq(1:seg(n).len) = f03_seq_est( seg(n).cvg_dep(:,1:seg(n).len) );
                            seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                        else
                            seg(n).len = seg(n).len - (tail_len_b - Min_tail_length);
                            tail_ind(n) = Min_tail_length;
                            seg(n).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(4:-1:1,(tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                            seg(n).seq(1:seg(n).len) = 3-seg(n).seq((tail_len_b - Min_tail_length)+seg(n).len:-1:(tail_len_b - Min_tail_length)+1);
                            seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ) );
                        end
                    end
                end
            end
        end
    end

    seg_cnt = 0;
    for n = 1:n_seg
        if b_seg_valid(n) == 1 && seg(n).len > num_N_mer
            seg_cnt = seg_cnt + 1;
            seg(seg_cnt).id = seg(n).id;
            seg(seg_cnt).len = seg(n).len;
            seg(seg_cnt).cvg_dep(:,1:seg(n).len) = seg(n).cvg_dep(:,1:seg(n).len);
            seg(seg_cnt).seq(1:seg(n).len) = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(n).len) );
            [axx, nc] = size(seg(n).cvg_dep(:,1:seg(n).len));
            if nc ~= seg(n).len || seg(n).len ~= length(seg(seg_cnt).seq(1:seg(n).len))
                str = sprintf('ERROR 3: %d ~= %d ~= %d', nc, seg(n).len, length(seg(seg_cnt).seq(1:seg(n).len)) );
                disp( str );
            end
            seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(n).len) ) );
            tail_ind(seg_cnt) = tail_ind(n);
        else
            b_seg_valid(n) = 0;
            cntg_valid_ind( seg(n).id ) = -2;
        end
    end
    n_seg = seg_cnt;
end
% end % for cphs       

%% Save data
% str_disp = sprintf('   Saving contig info. ' );
% fprintf( '\n%s', str_disp );
% fprintf( fp_log, '\n%s', str_disp );
% 
% % fname_out2 = sprintf('%s_js%s_c%d', out_file_prefix, tver, cfg.connection_threshold_short );
% fname_out2 = sprintf('%s', out_file_prefix );
% fname_txt = sprintf('%s.cntg2', fname_out2 );
% fp_t = fopen( fname_txt, 'wt' );
% 
% kmod = round(n_seg/8);
% for k = 1:1:n_seg
%     fprintf(fp_t,'>Contig2\t%d\t%d\t%d\t%d\n', k, seg(k).id, N_reads_valid, N_reads_total );
%     fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg(k).seq(1:seg(k).len) ) );
%     for m2 = 1:1:4
%         for m1 = 1:1:seg(k).len-1
%             fprintf(fp_t,'%d\t', seg(k).cvg_dep(m2,m1) );
%         end
%         m1 = seg(k).len;
%         fprintf(fp_t,'%d\n', seg(k).cvg_dep(m2,m1) );
%     end
%             
%     if mod(k,kmod) == 0
%         fprintf('.');
%     end
% end
% fclose(fp_t);
% 
% fname_out2 = sprintf('%s', out_file_prefix );
% fname_txt = sprintf('%s.rcmap2', fname_out2 );
% fp_t = fopen( fname_txt, 'wt' );
% fprintf(fp_t, '%d\n', cntg_cnt);
% kmod = round(cntg_cnt/6);
% for k = 1:1:cntg_cnt
%     fprintf(fp_t, '%d\n', cntg_valid_ind(k));
%     if mod(k,kmod) == 0
%         fprintf('.');
%     end
% end
% Nrt = round( N_reads_total/(R_mode+1) );
% fprintf(fp_t, '%d\n', Nrt);
% kmod = round(Nrt/6);
% for k = 1:1:Nrt
%     fprintf(fp_t, '%d\t%d\n', read_cntg_map(k,1), read_cntg_map(k,2) );
%     if mod(k,kmod) == 0
%         fprintf('.');
%     end
% end
% fclose(fp_t);
% 
% str_disp = sprintf(' saved to %s', fname_txt );
% fprintf( '%s ', str_disp );
% fprintf( fp_log, '%s ', str_disp );

end % soption

%% remove Segs having cvg depth less than the minimum value
if cfg.min_cvg_depth_js > 1
    str_disp = sprintf('   Removing contigs with cvg depth below %f ', cfg.min_cvg_depth_js );
    fprintf( '\n%s', str_disp );
    fprintf( fp_log, '\n%s', str_disp );

    seg_cnt = 0;
    n_step = round(n_seg/12);
    for n = 1:n_seg
        seg(n).seq = f03_seq_est( seg(n).cvg_dep(:,1:seg(n).len) );
        seg(n).ave_cvg_dep = mean( sum( seg(n).cvg_dep(:,1:seg(n).len) ));

        if seg(n).ave_cvg_dep > cfg.min_cvg_depth_js && seg(n).len > cfg.min_cntg_length_js 
            seg_cnt = seg_cnt + 1;
            if seg_cnt ~= n
                seg(seg_cnt).id = seg(n).id;
                seg(seg_cnt).len = seg(n).len;
                seg(seg_cnt).cvg_dep = seg(n).cvg_dep(:,1:seg(n).len);
                seg(seg_cnt).seq = f03_seq_est( seg(seg_cnt).cvg_dep(:,1:seg(n).len) );
                seg(seg_cnt).ave_cvg_dep = mean( sum( seg(seg_cnt).cvg_dep(:,1:seg(n).len) ) );
            end
        else
            cntg_valid_ind( seg(n).id ) = -2;
        end
        if mod(n, n_step) == 0
            fprintf('.');
        end
    end
    str_disp = sprintf('%d removed -> N_cntg: %d', n_seg-seg_cnt, seg_cnt );
    fprintf( '%s', str_disp );
    fprintf( fp_log, '%s', str_disp );
    n_seg = seg_cnt;
end

%% Junction search --> Create junction info.
fprintf( '\n   Searching junctions' );
fprintf( fp_log, '\n   Searching junctions' );
str_hdr = sprintf(' .....' );
Nchar = fprintf( '%s', str_hdr );
fprintf( fp_log, '%s', str_hdr );

cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth) cfg.connection_threshold_short cfg.bool_ss_ind];
num_N_mer = 8;
seg_sel_dsr = 4;
Max_overlap_window = cfg.safe_overlap_threshold; 
Num_proc_sim = cfg.proc_unit_size*2;
N_Loops = ceil( n_seg/ Num_proc_sim );
max_ovlp = 0;
min_ovlp = 1000;
num_divided = zeros( n_seg, 1 );
group_id = zeros(n_seg, 1);
N_groups = 0;
N_junctions = 0;
junction_info = int32( zeros(n_seg*4, 14) );
seg_connection_mode_lst = zeros( n_seg*100, 6, 'int32' );
seg_connection_mode_cnt = 0;
n_seg_sel = 0;

%% Update cntg_valid_ind
if R_mode > 0
    
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) < 0
        % cntg_gid(k) = 0;
    else
        if cntg_valid_ind(k) == 0 || cntg_valid_ind(k) == k
            cntg_valid_ind(k) = k;
        else
            cid_p = k;
            gid_tmp = 0;
            while(1)
                cid = cntg_valid_ind(cid_p);
                if cntg_valid_ind(cid) < 0
                    gid_tmp = -1;
                    break;
                else
                    if cntg_valid_ind(cid) == 0 || cntg_valid_ind(cid) == cid
                        gid_tmp = cid;
                        break;
                    else
                        cid_p = cid;
                    end
                end
            end
            cid_p = k;
            while(1)
                cid = cntg_valid_ind(cid_p);
                cntg_valid_ind(cid_p) = gid_tmp;
                if cntg_valid_ind(cid) <= 0 || cntg_valid_ind(cid) == cid
                    cntg_valid_ind(cid) = gid_tmp;
                    break;
                else
                    cid_p = cid;
                end
            end
        end
    end
end

%% set cid_to_seg_map using cntg_valid_ind
cid_to_seg_map = zeros( cntg_cnt, 1, 'uint32' );
for k = 1:1:n_seg
    cid_to_seg_map( seg(k).id ) = k;
end
N_error = 0;
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) > 0 % && cntg_valid_ind(k) ~= k
        if cid_to_seg_map( k ) == 0
            cid_to_seg_map( k ) = cid_to_seg_map( cntg_valid_ind(k) );
            if cntg_valid_ind( cntg_valid_ind(k) ) ~= cntg_valid_ind(k)
                N_error = N_error + 1;
            end
        else
            if cid_to_seg_map( k ) ~= cid_to_seg_map( cntg_valid_ind(k) )
                N_error = N_error + 1;
            end
        end
    end
end
if N_error > 0
    fprintf(' ERROR: set cid_to_seg_map .. %d ', N_error );
end

%% Set list of possible junctions - Use read_cntg_map, cntg_valid_ind, cid_to_seg_map
Max_n_junctions_per_cntg = 200;
jpc_tmp.bs = Max_n_junctions_per_cntg;
jpc_tmp.nj = 0;
jpc_tmp.ji = zeros( Max_n_junctions_per_cntg, 1, 'uint32' );
jpc = repmat( jpc_tmp, n_seg, 1 );
n_pj = 0;
% if R_mode > 0
Nrt = round( N_reads_total/(R_mode+1) );
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 && read_cntg_map(k,2) > 0
        if read_cntg_map(k,1) ~= read_cntg_map(k,2)
            cntg1 = cntg_valid_ind( read_cntg_map(k,1) );
            cntg2 = cntg_valid_ind( read_cntg_map(k,2) );
            if cntg1 == 0
                cntg1 = read_cntg_map(k,1);
            end
            if cntg2 == 0
                cntg2 = read_cntg_map(k,2);
            end
            if cntg1 > 0 && cntg2 > 0
                if cntg1 ~= cntg2
                    for m = 1:1:2
                        if m == 1
                            cid1 = cid_to_seg_map( cntg1 );
                            cid2 = cid_to_seg_map( cntg2 );
                        else
                            cid2 = cid_to_seg_map( cntg1 );
                            cid1 = cid_to_seg_map( cntg2 );
                        end
                        b_tmp = 1;
                        if jpc( cid1 ).nj > 0
                            jmatch = find( jpc( cid1 ).ji(1:jpc( cid1 ).nj) == cid2, 1 );
                            if ~isempty( jmatch )
                                b_tmp = 0;
                            end
                        end
                        if b_tmp > 0
                            jpc( cid1 ).nj = jpc( cid1 ).nj + 1;
                            if jpc( cid1 ).nj > jpc( cid1 ).bs
                                jpc( cid1 ).ji = [jpc( cid1 ).ji; zeros( Max_n_junctions_per_cntg, 1, 'uint32' ) ];
                                jpc( cid1 ).bs = jpc( cid1 ).bs + Max_n_junctions_per_cntg;
                            end
                            jpc( cid1 ).ji( jpc( cid1 ).nj ) = cid2;
                            if m == 1
                                n_pj = n_pj +1;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% pre-junction search

    cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth) cfg.connection_threshold_short cfg.bool_ss_ind];
    % s_idxs = find( [jpc(1:end).nj] > 0 );
    for n = 1:1:n_seg % length(s_idxs) % n_seg
        dividend_idx = n; %s_idxs(n); 

        if group_id(dividend_idx) == 0
            N_groups = N_groups + 1;
            group_id(dividend_idx) = N_groups;
        else
        end

        if jpc(dividend_idx).nj > 0 && seg(dividend_idx).len > cfg.short_seg_threshold %cfg.nominal_read_length

            % tcfg = cfg;
            % Overlap_threshold = cfg.connection_threshold;
            % Ave_num_error_normalized = ceil( cfg.norm_dist_threshold*(Overlap_threshold) );
            % Match_threshold = round( max( (Overlap_threshold - (num_N_mer-seg_sel_dsr))/seg_sel_dsr - Ave_num_error_normalized*(num_N_mer/seg_sel_dsr), Min_sel_th ) );

            if seg(dividend_idx).len >= cfg.short_seg_threshold 
                
                blk_len = Blk_len; 
                blk_len = blk_len - mod( blk_len, 2 );
                blk_not_ovlp = blk_len - round(cfg.nominal_read_length*1.1); 
                if seg( dividend_idx ).len <=  blk_len
                    n_blk = 1;
                else
                    n_blk = 0;
                    n_blk_max = ceil( seg( dividend_idx ).len/blk_not_ovlp );
                    for m = 1:1:n_blk_max
                        if blk_not_ovlp*(m-1)+blk_len >= seg( dividend_idx ).len
                            n_blk = m;
                            break
                        end
                    end
                    if n_blk == 0
                        fprintf('    ERROR: n_blk ');
                    end
                end

                n_js = 0;
                % b_cmplt = zeros(2,1);
                % [seg_seq_Nmer, seg_len_Nmer] = sub_NumSeq2NmerSeq( seg(dividend_idx).seq, num_N_mer );

                j_cnt_divisor = zeros(n_seg,2);
                j_idx_divisor = zeros(100,n_seg,2);

                for mn = 1:1:n_blk

                    if n_blk == 1
                        startp = 1;
                        endp = seg(dividend_idx).len;
                    else
                        startp = blk_not_ovlp*(mn-1)+1;
                        if mn == n_blk
                            endp = seg(dividend_idx).len;
                        else
                            endp = blk_not_ovlp*(mn-1)+blk_len;
                        end
                    end
                    % seg_seq_Nmer_tmp = seg_seq_Nmer( startp:endp-num_N_mer+1 );
                    
                    for kht = 1:1:2

%                         if kht == 1
%                             Nmer_Matches = sum( Seg_Nmer_map_head( seg_idx_valid, seg_seq_Nmer_tmp(1:seg_sel_dsr:end)+1 ), 2 );
%                         else
%                             Nmer_Matches = sum( Seg_Nmer_map_tail( seg_idx_valid, seg_seq_Nmer_tmp(1:seg_sel_dsr:end)+1 ), 2 );
%                         end
%                         siv = seg_idx_valid( Nmer_Matches >= Match_threshold );
%                         if ~isempty(siv)
%                             [sorted_matches, ssi] = sort( Nmer_Matches(siv), 'descend' );
%                             sel_seg_idx = siv(ssi);
%                         else
%                             sel_seg_idx = siv;
%                         end
                        sel_seg_idx = jpc(dividend_idx).ji(1:jpc(dividend_idx).nj);
                        sp = 1;
                        n_seg_sel = 0.99*n_seg_sel + 0.01*length(sel_seg_idx);

                        for k = 1:1:length(sel_seg_idx) %min( Max_num_sel, length(sel_seg_idx) ) 
                            divisor_idx = sp + sel_seg_idx(k) -1; 

                            if seg( divisor_idx ).len >= Min_seg_length && divisor_idx ~= dividend_idx 

                                junc_mode = [0 0];
                                junction_info_tmp = int32( zeros(2, 14) );
                                % if b_cmplt(1) == 0
%                                 fprintf('(%d,%d)\n', dividend_idx, divisor_idx ); %, ...
%                                 fprintf('(%d,%d) - %s, %s\n', dividend_idx, divisor_idx, ...
%                                     sub_NumSeq2NTStr( seg(dividend_idx).seq(startp:endp) ), sub_NumSeq2NTStr( seg(divisor_idx).seq ) );
                                    [ovlp_pos, ovlp_len, j_mode, dst, njn] = sub_tjn( ...
                                        seg(dividend_idx).seq(startp:endp), seg(divisor_idx).seq, cfg_val, kht-1, Run_thresh );
                                    % junction_mode = 1 -> divisor head, split
                                    % junction_mode = 2 -> divisor tail, merge 
                                    % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                    % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 

                                    c_tmp = 0;
                                    if (mn < n_blk)
                                        if (ovlp_pos > blk_not_ovlp)
                                            c_tmp = 1;
                                        end
                                    end
                                    if (mn > 1) && (ovlp_pos == 1)
                                        c_tmp = 1;
                                    end
                                    if (ovlp_pos ~= 0) && c_tmp ==0
                                        if ovlp_pos < 0
                                            if b_inter_disp > 0
                                                fprintf('\n      ERROR JD: error occurred, %d ', njn);
                                            end
                                        end

                                        % slice dividend
                                        ovlp_pos = ovlp_pos + startp-1;
                                        bndry1 = ovlp_pos; 
                                        bndry2 = ovlp_pos + ovlp_len -1;
                                        if max_ovlp < ovlp_len
                                            max_ovlp = ovlp_len;
                                        end
                                        if min_ovlp > ovlp_len
                                            min_ovlp = ovlp_len;
                                        end
                                        % junction_mode = 1 -> divisor head, split
                                        % junction_mode = 2 -> divisor tail, merge 
                                        % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                        % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
                                        if bndry1 > 1 && bndry2 < seg(dividend_idx).len
                                            Jnob = Junction_overlap_backoff; 
                                            switch( j_mode )
                                                case 1,
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                case 2,
                                                    ovlp_pos = ovlp_pos + Jnob;
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                case 3,
                                                    ovlp_pos = ovlp_pos + Jnob;
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                otherwise,
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                            end
                                        end
                                        b_tmp = 0;
                                        if bndry1 < bndry2
                                            if bndry1 == 1 || bndry2 == seg(dividend_idx).len % if this is connection
                                                if dividend_idx < divisor_idx
                                                    b_tmp = 1;
                                                else
                                                end
                                            end
                                        else
                                            b_tmp = 1;
                                        end

                                        if b_tmp == 0
                                            % junction_mode = 1 -> divisor head, split
                                            % junction_mode = 2 -> divisor tail, merge 
                                            % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                            % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
                                            junction_info_tmp(kht,1) = group_id(dividend_idx);
                                            junction_info_tmp(kht,2) = dividend_idx;
                                            junction_info_tmp(kht,3) = divisor_idx;

                                            junction_info_tmp(kht,5) = ovlp_len;
                                            junction_info_tmp(kht,6) = ovlp_pos;
                                            junction_info_tmp(kht,7) = ovlp_pos + ovlp_len -1;
                                            junction_info_tmp(kht,8) = seg(dividend_idx).len;

                                            junction_info_tmp(kht,11) = seg(divisor_idx).len;
                                            junction_info_tmp(kht,12) = j_mode;
                                            junction_info_tmp(kht,14) = dst;

                                            junc_mode(kht) = j_mode;
                                            switch( j_mode )
                                                case 1,
                                                    junction_info_tmp(kht,4) = 1; 
                                                    junction_info_tmp(kht,9) = 1;
                                                    junction_info_tmp(kht,10) = ovlp_len;
                                                    if kht == 2 && b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 03');
                                                    end
                                                case 2,
                                                    junction_info_tmp(kht,4) = 1; 
                                                    junction_info_tmp(kht,9) = seg(divisor_idx).len - ovlp_len + 1;
                                                    junction_info_tmp(kht,10) = seg(divisor_idx).len;
                                                    if kht == 1 && b_inter_disp > 0
                                                        fprintf('\n      ERROR JD: This case cannot happen 01');
                                                    end
                                                case 3,
                                                    junction_info_tmp(kht,4) = -1; 
                                                    junction_info_tmp(kht,9) = 1; 
                                                    junction_info_tmp(kht,10) = ovlp_len; 
                                                    if kht == 2 b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 04');
                                                    end
                                                otherwise % case 4
                                                    junction_info_tmp(kht,4) = -1; 
                                                    junction_info_tmp(kht,9) = seg(divisor_idx).len - ovlp_len + 1;
                                                    junction_info_tmp(kht,10) = seg(divisor_idx).len;
                                                    if kht == 1 && b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 02');
                                                    end
                                            end
                                            num_divided( dividend_idx ) = num_divided( dividend_idx ) + 1;
                                        end
                                    end
                                % end

                                if junc_mode(1) ~= 0 || junc_mode(2) ~= 0
                                    % check compatibility
                                    [b_compatibility, cm_lst] = f04_compatibility_check( junction_info_tmp(1,:), junction_info_tmp(2,:), b_inter_disp );
                                    if b_compatibility(1) == 0 && b_compatibility(2) == 0
                                    else
                                        % set group id
                                        if group_id(divisor_idx) == 0 
                                            group_id(divisor_idx) = group_id(dividend_idx);
                                        else
                                            if group_id(dividend_idx) == group_id(divisor_idx)
                                            else
                                                gid_tmp = group_id(divisor_idx);
                                                jidxs = find( junction_info(1:N_junctions,1) == gid_tmp );
                                                if ~isempty(jidxs)
                                                    junction_info(jidxs,1) = group_id(dividend_idx);
                                                end
                                                sidxs = find( group_id(1:n_seg) == gid_tmp );
                                                if ~isempty(sidxs)
                                                    group_id(sidxs) = group_id(dividend_idx);
                                                end
                                            end
                                        end
                                        % add to junction_info_cur
                                        for m1 = 1:1:2
                                            if junc_mode(m1) ~= 0 && b_compatibility(m1) ~= 0
                                                junction_info_tmp(m1,1) = group_id(dividend_idx);
                                                junction_info_tmp(m1,13) = b_compatibility(m1);
                                                n_js = n_js + 1;
                                                N_junctions = N_junctions +1;
                                                junction_info(N_junctions,:) = junction_info_tmp(m1,:);
                                                j_cnt_divisor(divisor_idx,m1) = j_cnt_divisor(divisor_idx,m1) + 1;
                                                j_idx_divisor(j_cnt_divisor(divisor_idx,m1),divisor_idx,m1) = N_junctions;
                                            end
                                        end
                                    end
                                end
                                % if sum( b_cmplt ) == 2
                                %     break;
                                % end
                           end 
                        end % for k
                    end % kht
                end % mn
                n_possible_incompatible_junctions = 0;
                s_idx_t = find(sum( j_cnt_divisor, 2 ));
                for k = 1:1:length(s_idx_t)
                    divisor_idx = s_idx_t(k); 
                    b_m_fnd = j_cnt_divisor(divisor_idx,:)';
                   if b_m_fnd(1) > 0 || b_m_fnd(2) > 0 
                        % Select only one
                        if b_m_fnd(1) > 0 && b_m_fnd(2) > 0
                            % Test overlap length
                            ovlp2 = zeros(b_m_fnd(1),b_m_fnd(2));
                            for m1 = 1:1:b_m_fnd(1)
                                ji_1 = j_idx_divisor(m1,divisor_idx,1);
                                for m2 = 1:1:b_m_fnd(2)
                                    ji_2 = j_idx_divisor(m2,divisor_idx,2);
                                    [b_comp, axx] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                    if sum(b_comp) > 0
                                        ov1 = junction_info(ji_1,5)-junction_info(ji_1,14);
                                        ov2 = junction_info(ji_2,5)-junction_info(ji_2,14);
                                        ovlp2(m1,m2) = ov1*b_comp(1) + ov2*b_comp(2);
                                    end
                                end
                            end
                            if b_m_fnd(1) == 1
                                mri(1) = 1;
                                [mxv(1),mci(1)] = max( ovlp2 );
                                % mxv(1) = ovlp2( mri(1), mci(1) );
                                ovlp2( mri(1), mci(1) ) = 0;
                                mri(2) = 1;
                                [mxv(2),mci(2)] = max( ovlp2 );
                                % mxv(2) = ovlp2( mri(2), mci(2) );
                            else
                                [mxvs, mris] = max(ovlp2);
                                [mxv(1),mci(1)] = max( mxvs );
                                mri(1) = mris(mci(1));
                                % mxv(1) = ovlp2( mri(1), mci(1) );
                                ovlp2( mri(1), mci(1) ) = 0;
                                [mxvs2, mris2] = max(ovlp2);
                                [mxv(2),mci(2)] = max( mxvs2 );
                                mri(2) = mris2(mci(2));
                                % mxv(2) = ovlp2( mri(2), mci(2) );
                            end
                            if mxv(1) > mxv(2) || mxv(2) == 0
                                ji_1 = j_idx_divisor(mri(1),divisor_idx,1);
                                ji_2 = j_idx_divisor(mci(1),divisor_idx,2);
                                if ji_1 == 0 || ji_2 == 0
                                    fprintf('\n      ERROR: ji_1 == 0 || ji_2 == 0 (A) for (%d,%d)', divisor_idx, dividend_idx );
                                else
                                    [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                    if b_comp(1) > 0
                                        for mx = 1:1:b_m_fnd(1)
                                            if mx ~= mri(1)
                                                ji = j_idx_divisor(mx,divisor_idx,1);
                                                junction_info(ji,13) = 0;
                                            end
                                        end
                                    end
                                    if b_comp(2) > 0
                                        for mx = 1:1:b_m_fnd(2)
                                            if mx ~= mci(1)
                                                ji = j_idx_divisor(mx,divisor_idx,2);
                                                junction_info(ji,13) = 0;
                                            end
                                        end
                                    end
                                    % if b_comp(1) > 0
                                    %     N_junctions = N_junctions +1;
                                    %     junction_info(N_junctions,:) = junction_info_cur(:,mri(1),1)';
                                    % end
                                    % if b_comp(2) > 0
                                    %     N_junctions = N_junctions +1;
                                    %     junction_info(N_junctions,:) = junction_info_cur(:,mci(1),2)';
                                    % end
                                    if isempty( cm_lst )
                                    else
                                        if junction_info(ji_1,5)-junction_info(ji_1,14) > junction_info(ji_2,5)-junction_info(ji_2,14)
                                            boff = ji_1;
                                        else
                                            boff = ji_2;
                                        end
                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) boff];
                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) boff];
                                    end
                                    b_m_fnd(1) = 1;
                                    b_m_fnd(2) = 1;
                                end
                            else
                                for mz = 1:1:1
                                    ji_1 = j_idx_divisor(mri(mz),divisor_idx,1);
                                    ji_2 = j_idx_divisor(mci(mz),divisor_idx,2);
                                    if ji_1 == 0 || ji_2 == 0
                                        fprintf('\n      ERROR: ji_1 == 0 || ji_2 == 0 (B) for (%d,%d)', divisor_idx, dividend_idx );
                                    else
                                        [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                        % leave them
                                        if b_comp(1) > 0
                                            for mx = 1:1:b_m_fnd(1)
                                                if mx ~= mri(mz)
                                                    ji = j_idx_divisor(mx,divisor_idx,1);
                                                    junction_info(ji,13) = 0;
                                                end
                                            end
                                            b_m_fnd(1) = 1;
                                        end
                                        if b_comp(2) > 0
                                            for mx = 1:1:b_m_fnd(2)
                                                if mx ~= mri(mz)
                                                    ji = j_idx_divisor(mx,divisor_idx,2);
                                                    junction_info(ji,13) = 0;
                                                end
                                            end
                                            b_m_fnd(2) = 1;
                                        end
                                        if isempty( cm_lst )
                                        else
                                            if junction_info(ji_1,5)-junction_info(ji_1,14) > junction_info(ji_2,5)-junction_info(ji_2,14)
                                                boff = ji_1;
                                            else
                                                boff = ji_2;
                                            end
                                            seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                            seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) boff];
                                            seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                            seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) boff];
                                        end
                                    end
                                end
                            end
                        else % b_m_fnd(1) == 0 || b_m_fnd(2) == 0
                            for mz = 1:1:2
                                if b_m_fnd(mz) > 0 
                                    ovlp = zeros(b_m_fnd(mz),1);
                                    for m1 = 1:1:b_m_fnd(mz)
                                        jit = j_idx_divisor(m1,divisor_idx,mz);
                                        ovlp(m1) = junction_info(jit,5)-junction_info(jit,14);
                                    end
                                    [mxv(1), mxi(1)] = max(ovlp);
                                    ovlp(mxi(1)) = 0;
                                    [mxv(2), mxi(2)] = max(ovlp);
                                    if mxv(1) > mxv(2) || mxv(2) == 0
                                        ji = j_idx_divisor(mxi(1),divisor_idx,mz);
                                        ji_tmp = zeros(1,14);
                                        if mz == 1
                                            [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji,:), ji_tmp, b_inter_disp );
                                        else
                                            [b_comp, cm_lst] = f04_compatibility_check( ji_tmp, junction_info(ji,:), b_inter_disp );
                                        end
                                        if sum( b_comp ) > 0
                                            for mx = 1:1:b_m_fnd(mz)
                                                if mx ~= mxi(1)
                                                    jit = j_idx_divisor(mx,divisor_idx,mz);
                                                    junction_info(jit,13) = 0;
                                                end
                                            end
                                            if isempty( cm_lst )
                                            else
                                                seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                    [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) ji];
                                                seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                    [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) ji];
                                            end
                                            b_m_fnd(mz) = 1;
                                        end
                                    else
                                        ji_tmp = zeros(1,14);
                                        if mz == 1
                                            if mxi(1) == mxi(2)
                                                b_m_fnd(mz) = 1;
                                            else
                                                b_m_fnd(mz) = 1; 
                                            end
                                            for my = 1:1:b_m_fnd(mz)
                                                ji = j_idx_divisor(mxi(my),divisor_idx,mz);
                                                [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji,:), ji_tmp, b_inter_disp );
                                                if sum( b_comp ) > 0
                                                    for mx = 1:1:b_m_fnd(mz)
                                                        if mx ~= mxi(my)
                                                            jit = j_idx_divisor(mx,divisor_idx,mz);
                                                            junction_info(jit,13) = 0;
                                                        end
                                                    end
                                                    if isempty( cm_lst )
                                                    else
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(my) ji];
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(my) ji];
                                                    end
                                                else
                                                end
                                            end
                                        else % if mz == 2
                                            if mxi(1) == mxi(2)
                                                b_m_fnd(mz) = 1;
                                            else
                                                b_m_fnd(mz) = 1; 
                                            end
                                            for my = 1:1:b_m_fnd(mz)
                                                ji = j_idx_divisor(mxi(my),divisor_idx,mz);
                                                [b_comp, cm_lst] = f04_compatibility_check( ji_tmp, junction_info(ji,:), b_inter_disp );
                                                if sum( b_comp ) > 0
                                                    for mx = 1:1:b_m_fnd(mz)
                                                        if mx ~= mxi(my)
                                                            jit = j_idx_divisor(mx,divisor_idx,mz);
                                                            junction_info(jit,13) = 0;
                                                        end
                                                    end
                                                    if isempty( cm_lst )
                                                    else
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(my) ji]; 
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(my) ji]; 
                                                    end
                                                end
                                            end
                                        end
                                    end
                                else
                                    b_m_fnd(mz) = 0;
                                end
                            end
                        end
                    end % if b_m_fnd(1) > 1 || b_m_fnd(2) > 1 %n_js > 2
                    if b_m_fnd(1) > 1 || b_m_fnd(2) > 1 %n_js > 2
                        if b_inter_disp > 0
                            str_disp = sprintf('      WARNING: (%d,%d) junctions for (%d, %d)  ', b_m_fnd(1), b_m_fnd(2), dividend_idx, divisor_idx );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            Nchar = 0;
                        end
                    end
                    n_possible_incompatible_junctions = n_possible_incompatible_junctions + max( b_m_fnd(1)-1, 0 ) + max( b_m_fnd(2)-1, 0 );        
                end
                % if mod(n, n_step) == 0
                    % fprintf('.');
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = fprintf( '%s %d/%d/%d ', str_hdr, n_seg, n, N_junctions ); 
                % end
            end
        end % if seg(dividend_idx).len    
    end % for n = 1:1:n_seg
else

%% Junction search start
Nchar = 0;
% if Nchar > 0
cfg_val = [cfg.norm_dist_threshold cfg.connection_threshold_short round(cfg.max_overlap_depth) cfg.connection_threshold cfg.bool_ss_ind];
for lp = 1:1:N_Loops
    
    sp = (lp-1)*Num_proc_sim+1;
    if lp < N_Loops
        ep = lp*Num_proc_sim;
    else
        ep = n_seg;
    end
    num_proc_sim = ep-sp+1;

    b_seg_length = zeros(num_proc_sim,1);
    Seg_Nmer_map_head = zeros( num_proc_sim, 4^num_N_mer, 'int8' );
    Seg_Nmer_map_tail = zeros( num_proc_sim, 4^num_N_mer, 'int8' );
    for n = 1:num_proc_sim 
        idx = sp+n-1; 
        SLen = seg(idx).len;
        if SLen > Max_overlap_window 
            b_seg_length(n) = 1;
            OLen = min( Max_overlap_window, seg(idx).len );
            [Nmer1, sl1] = sub_NumSeq2NmerSeq( seg(idx).seq(1:OLen), num_N_mer );
            [Nmer2, sl2] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(OLen:-1:1), num_N_mer );
            [Nmer3, sl3] = sub_NumSeq2NmerSeq( seg(idx).seq(SLen-OLen+1:SLen), num_N_mer );
            [Nmer4, sl4] = sub_NumSeq2NmerSeq( 3-seg(idx).seq(SLen:-1:SLen-OLen+1), num_N_mer );
            Seg_Nmer_map_head( n, [Nmer1(1:sl1) Nmer2(1:sl2)]+1 ) = 1;
            Seg_Nmer_map_tail( n, [Nmer3(1:sl3) Nmer4(1:sl4)]+1 ) = 1;
        end
        % if mod(n, n_step) == 0
        %     fprintf('.' );
        %     % Nchar = Nchar + 1; % fprintf( '%d/%d', (n_seg - n), n_seg ); 
        % end
    end
    seg_idx_valid = find( b_seg_length );
    % Seg_Nmer_map_tmp = int8( zeros( 1, 4^num_N_mer ) );

    for n = 1:1:n_seg
        dividend_idx = n; 

        if group_id(dividend_idx) == 0
            N_groups = N_groups + 1;
            group_id(dividend_idx) = N_groups;
        else
        end

        if seg(dividend_idx).len > cfg.short_seg_threshold %cfg.nominal_read_length

            % tcfg = cfg;
            Overlap_threshold = cfg.connection_threshold;
            Ave_num_error_normalized = ceil( cfg.norm_dist_threshold*(Overlap_threshold) );
            Match_threshold = round( max( (Overlap_threshold - (num_N_mer-seg_sel_dsr))/seg_sel_dsr - Ave_num_error_normalized*(num_N_mer/seg_sel_dsr), Min_sel_th ) );

            if seg(dividend_idx).len >= cfg.short_seg_threshold 
                
                blk_len = Blk_len; 
                blk_len = blk_len - mod( blk_len, 2 );
                blk_not_ovlp = blk_len - round(cfg.nominal_read_length*1.1); 
                if seg( dividend_idx ).len <=  blk_len
                    n_blk = 1;
                else
                    n_blk = 0;
                    n_blk_max = ceil( seg( dividend_idx ).len/blk_not_ovlp );
                    for m = 1:1:n_blk_max
                        if blk_not_ovlp*(m-1)+blk_len >= seg( dividend_idx ).len
                            n_blk = m;
                            break
                        end
                    end
                    if n_blk == 0
                        fprintf('    ERROR: n_blk ');
                    end
                end

                n_js = 0;
                % b_cmplt = zeros(2,1);
                [seg_seq_Nmer, seg_len_Nmer] = sub_NumSeq2NmerSeq( seg(dividend_idx).seq, num_N_mer );

                j_cnt_divisor = zeros(n_seg,2);
                j_idx_divisor = zeros(100,n_seg,2);

                for mn = 1:1:n_blk

                    if n_blk == 1
                        startp = 1;
                        endp = seg(dividend_idx).len;
                    else
                        startp = blk_not_ovlp*(mn-1)+1;
                        if mn == n_blk
                            endp = seg(dividend_idx).len;
                        else
                            endp = blk_not_ovlp*(mn-1)+blk_len;
                        end
                    end
                    seg_seq_Nmer_tmp = seg_seq_Nmer( startp:endp-num_N_mer+1 );
                    
                    for kht = 1:1:2

                        if kht == 1
                            Nmer_Matches = sum( Seg_Nmer_map_head( seg_idx_valid, seg_seq_Nmer_tmp(1:seg_sel_dsr:end)+1 ), 2 );
                        else
                            Nmer_Matches = sum( Seg_Nmer_map_tail( seg_idx_valid, seg_seq_Nmer_tmp(1:seg_sel_dsr:end)+1 ), 2 );
                        end
%                         siv = seg_idx_valid( Nmer_Matches >= Match_threshold );
%                         if ~isempty(siv)
%                             [sorted_matches, ssi] = sort( Nmer_Matches(siv), 'descend' );
%                             sel_seg_idx = siv(ssi);
%                         else
%                             sel_seg_idx = siv;
%                         end
                        siv = find( Nmer_Matches >= Match_threshold );
                        if ~isempty(siv)
                            [sorted_matches, ssi] = sort( Nmer_Matches(siv), 'descend' );
                            sel_seg_idx = seg_idx_valid(siv(ssi));
                        else
                            sel_seg_idx = siv;
                        end
                        
                        n_seg_sel = 0.99*n_seg_sel + 0.01*length(sel_seg_idx);

                        for k = 1:1:min( Max_num_sel, length(sel_seg_idx) ) 
                            divisor_idx = sp + sel_seg_idx(k) -1; 

                            if seg( divisor_idx ).len >= Min_seg_length && divisor_idx ~= dividend_idx 

                                junc_mode = [0 0];
                                junction_info_tmp = int32( zeros(2, 14) );
                                % if b_cmplt(1) == 0
                                    [ovlp_pos, ovlp_len, j_mode, dst, njn] = sub_tjn( ...
                                        seg(dividend_idx).seq(startp:endp), seg(divisor_idx).seq, cfg_val, kht-1, Run_thresh );
                                    % junction_mode = 1 -> divisor head, split
                                    % junction_mode = 2 -> divisor tail, merge 
                                    % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                    % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 

                                    c_tmp = 0;
                                    if (mn < n_blk)
                                        if (ovlp_pos > blk_not_ovlp)
                                            c_tmp = 1;
                                        end
                                    end
                                    if (mn > 1) && (ovlp_pos == 1)
                                        c_tmp = 1;
                                    end
                                    if (ovlp_pos ~= 0) && c_tmp ==0
                                        if ovlp_pos < 0
                                            if b_inter_disp > 0
                                                fprintf('\n      ERROR JD: error occurred, %d ', njn);
                                            end
                                        end

                                        % slice dividend
                                        ovlp_pos = ovlp_pos + startp-1;
                                        bndry1 = ovlp_pos; 
                                        bndry2 = ovlp_pos + ovlp_len -1;
                                        if max_ovlp < ovlp_len
                                            max_ovlp = ovlp_len;
                                        end
                                        if min_ovlp > ovlp_len
                                            min_ovlp = ovlp_len;
                                        end
                                        % junction_mode = 1 -> divisor head, split
                                        % junction_mode = 2 -> divisor tail, merge 
                                        % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                        % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
                                        if bndry1 > 1 && bndry2 < seg(dividend_idx).len
                                            Jnob = Junction_overlap_backoff; 
                                            switch( j_mode )
                                                case 1,
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                case 2,
                                                    ovlp_pos = ovlp_pos + Jnob;
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                case 3,
                                                    ovlp_pos = ovlp_pos + Jnob;
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                                otherwise,
                                                    ovlp_len = ovlp_len - Jnob;
                                                    bndry1 = ovlp_pos; 
                                                    bndry2 = ovlp_pos + ovlp_len -1;
                                            end
                                        end
                                        b_tmp = 0;
                                        if bndry1 < bndry2
                                            if bndry1 == 1 || bndry2 == seg(dividend_idx).len % if this is connection
                                                if dividend_idx < divisor_idx
                                                    b_tmp = 1;
                                                else
                                                end
                                            end
                                        else
                                            b_tmp = 1;
                                        end

                                        if b_tmp == 0
                                            % junction_mode = 1 -> divisor head, split
                                            % junction_mode = 2 -> divisor tail, merge 
                                            % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
                                            % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
                                            junction_info_tmp(kht,1) = group_id(dividend_idx);
                                            junction_info_tmp(kht,2) = dividend_idx;
                                            junction_info_tmp(kht,3) = divisor_idx;

                                            junction_info_tmp(kht,5) = ovlp_len;
                                            junction_info_tmp(kht,6) = ovlp_pos;
                                            junction_info_tmp(kht,7) = ovlp_pos + ovlp_len -1;
                                            junction_info_tmp(kht,8) = seg(dividend_idx).len;

                                            junction_info_tmp(kht,11) = seg(divisor_idx).len;
                                            junction_info_tmp(kht,12) = j_mode;
                                            junction_info_tmp(kht,14) = dst;

                                            junc_mode(kht) = j_mode;
                                            switch( j_mode )
                                                case 1,
                                                    junction_info_tmp(kht,4) = 1; 
                                                    junction_info_tmp(kht,9) = 1;
                                                    junction_info_tmp(kht,10) = ovlp_len;
                                                    if kht == 2 && b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 03');
                                                    end
                                                case 2,
                                                    junction_info_tmp(kht,4) = 1; 
                                                    junction_info_tmp(kht,9) = seg(divisor_idx).len - ovlp_len + 1;
                                                    junction_info_tmp(kht,10) = seg(divisor_idx).len;
                                                    if kht == 1 && b_inter_disp > 0
                                                        fprintf('\n      ERROR JD: This case cannot happen 01');
                                                    end
                                                case 3,
                                                    junction_info_tmp(kht,4) = -1; 
                                                    junction_info_tmp(kht,9) = 1; 
                                                    junction_info_tmp(kht,10) = ovlp_len; 
                                                    if kht == 2 b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 04');
                                                    end
                                                otherwise % case 4
                                                    junction_info_tmp(kht,4) = -1; 
                                                    junction_info_tmp(kht,9) = seg(divisor_idx).len - ovlp_len + 1;
                                                    junction_info_tmp(kht,10) = seg(divisor_idx).len;
                                                    if kht == 1 && b_inter_disp > 0
                                                        fprintf('\n      ERROR: This case cannot happen 02');
                                                    end
                                            end
                                            num_divided( dividend_idx ) = num_divided( dividend_idx ) + 1;
                                        end
                                    end
                                % end

                                if junc_mode(1) ~= 0 || junc_mode(2) ~= 0
                                    % check compatibility
                                    [b_compatibility, cm_lst] = f04_compatibility_check( junction_info_tmp(1,:), junction_info_tmp(2,:), b_inter_disp );
                                    if b_compatibility(1) == 0 && b_compatibility(2) == 0
                                    else
                                        % set group id
                                        if group_id(divisor_idx) == 0 
                                            group_id(divisor_idx) = group_id(dividend_idx);
                                        else
                                            if group_id(dividend_idx) == group_id(divisor_idx)
                                            else
                                                gid_tmp = group_id(divisor_idx);
                                                jidxs = find( junction_info(1:N_junctions,1) == gid_tmp );
                                                if ~isempty(jidxs)
                                                    junction_info(jidxs,1) = group_id(dividend_idx);
                                                end
                                                sidxs = find( group_id(1:n_seg) == gid_tmp );
                                                if ~isempty(sidxs)
                                                    group_id(sidxs) = group_id(dividend_idx);
                                                end
                                            end
                                        end
                                        % add to junction_info_cur
                                        for m1 = 1:1:2
                                            if junc_mode(m1) ~= 0 && b_compatibility(m1) ~= 0
                                                junction_info_tmp(m1,1) = group_id(dividend_idx);
                                                junction_info_tmp(m1,13) = b_compatibility(m1);
                                                n_js = n_js + 1;
                                                N_junctions = N_junctions +1;
                                                junction_info(N_junctions,:) = junction_info_tmp(m1,:);
                                                j_cnt_divisor(divisor_idx,m1) = j_cnt_divisor(divisor_idx,m1) + 1;
                                                j_idx_divisor(j_cnt_divisor(divisor_idx,m1),divisor_idx,m1) = N_junctions;
                                            end
                                        end
                                    end
                                end
                                % if sum( b_cmplt ) == 2
                                %     break;
                                % end
                           end 
                        end % for k
                    end % kht
                end % mn
                
                n_possible_incompatible_junctions = 0;
                s_idx_t = find(sum( j_cnt_divisor, 2 ));
                for k = 1:1:length(s_idx_t)
                    divisor_idx = s_idx_t(k); 
                    b_m_fnd = j_cnt_divisor(divisor_idx,:)';
                   if b_m_fnd(1) > 0 || b_m_fnd(2) > 0 
                        % Select only one
                        if b_m_fnd(1) > 0 && b_m_fnd(2) > 0
                            % Test overlap length
                            ovlp2 = zeros(b_m_fnd(1),b_m_fnd(2));
                            for m1 = 1:1:b_m_fnd(1)
                                ji_1 = j_idx_divisor(m1,divisor_idx,1);
                                for m2 = 1:1:b_m_fnd(2)
                                    ji_2 = j_idx_divisor(m2,divisor_idx,2);
                                    [b_comp, axx] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                    if sum(b_comp) > 0
                                        ov1 = junction_info(ji_1,5)-junction_info(ji_1,14);
                                        ov2 = junction_info(ji_2,5)-junction_info(ji_2,14);
                                        ovlp2(m1,m2) = ov1*b_comp(1) + ov2*b_comp(2);
                                    end
                                end
                            end
                            if b_m_fnd(1) == 1
                                mri(1) = 1;
                                [mxv(1),mci(1)] = max( ovlp2 );
                                % mxv(1) = ovlp2( mri(1), mci(1) );
                                ovlp2( mri(1), mci(1) ) = 0;
                                mri(2) = 1;
                                [mxv(2),mci(2)] = max( ovlp2 );
                                % mxv(2) = ovlp2( mri(2), mci(2) );
                            else
                                [mxvs, mris] = max(ovlp2);
                                [mxv(1),mci(1)] = max( mxvs );
                                mri(1) = mris(mci(1));
                                % mxv(1) = ovlp2( mri(1), mci(1) );
                                ovlp2( mri(1), mci(1) ) = 0;
                                [mxvs2, mris2] = max(ovlp2);
                                [mxv(2),mci(2)] = max( mxvs2 );
                                mri(2) = mris2(mci(2));
                                % mxv(2) = ovlp2( mri(2), mci(2) );
                            end
                            if mxv(1) > mxv(2) || mxv(2) == 0
                                ji_1 = j_idx_divisor(mri(1),divisor_idx,1);
                                ji_2 = j_idx_divisor(mci(1),divisor_idx,2);
                                if ji_1 == 0 || ji_2 == 0
                                    fprintf('\n      ERROR: ji_1 == 0 || ji_2 == 0 (A) for (%d,%d)', divisor_idx, dividend_idx );
                                else
                                    [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                    if b_comp(1) > 0
                                        for mx = 1:1:b_m_fnd(1)
                                            if mx ~= mri(1)
                                                ji = j_idx_divisor(mx,divisor_idx,1);
                                                junction_info(ji,13) = 0;
                                            end
                                        end
                                    end
                                    if b_comp(2) > 0
                                        for mx = 1:1:b_m_fnd(2)
                                            if mx ~= mci(1)
                                                ji = j_idx_divisor(mx,divisor_idx,2);
                                                junction_info(ji,13) = 0;
                                            end
                                        end
                                    end
                                    % if b_comp(1) > 0
                                    %     N_junctions = N_junctions +1;
                                    %     junction_info(N_junctions,:) = junction_info_cur(:,mri(1),1)';
                                    % end
                                    % if b_comp(2) > 0
                                    %     N_junctions = N_junctions +1;
                                    %     junction_info(N_junctions,:) = junction_info_cur(:,mci(1),2)';
                                    % end
                                    if isempty( cm_lst )
                                    else
                                        if junction_info(ji_1,5)-junction_info(ji_1,14) > junction_info(ji_2,5)-junction_info(ji_2,14)
                                            boff = ji_1;
                                        else
                                            boff = ji_2;
                                        end
                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) boff];
                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) boff];
                                    end
                                    b_m_fnd(1) = 1;
                                    b_m_fnd(2) = 1;
                                end
                            else
                                for mz = 1:1:1
                                    ji_1 = j_idx_divisor(mri(mz),divisor_idx,1);
                                    ji_2 = j_idx_divisor(mci(mz),divisor_idx,2);
                                    if ji_1 == 0 || ji_2 == 0
                                        fprintf('\n      ERROR: ji_1 == 0 || ji_2 == 0 (B) for (%d,%d)', divisor_idx, dividend_idx );
                                    else
                                        [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji_1,:), junction_info(ji_2,:), b_inter_disp );
                                        % leave them
                                        if b_comp(1) > 0
                                            for mx = 1:1:b_m_fnd(1)
                                                if mx ~= mri(mz)
                                                    ji = j_idx_divisor(mx,divisor_idx,1);
                                                    junction_info(ji,13) = 0;
                                                end
                                            end
                                            b_m_fnd(1) = 1;
                                        end
                                        if b_comp(2) > 0
                                            for mx = 1:1:b_m_fnd(2)
                                                if mx ~= mri(mz)
                                                    ji = j_idx_divisor(mx,divisor_idx,2);
                                                    junction_info(ji,13) = 0;
                                                end
                                            end
                                            b_m_fnd(2) = 1;
                                        end
                                        if isempty( cm_lst )
                                        else
                                            if junction_info(ji_1,5)-junction_info(ji_1,14) > junction_info(ji_2,5)-junction_info(ji_2,14)
                                                boff = ji_1;
                                            else
                                                boff = ji_2;
                                            end
                                            seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                            seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) boff];
                                            seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                            seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) boff];
                                        end
                                    end
                                end
                            end
                        else % b_m_fnd(1) == 0 || b_m_fnd(2) == 0
                            for mz = 1:1:2
                                if b_m_fnd(mz) > 0 
                                    ovlp = zeros(b_m_fnd(mz),1);
                                    for m1 = 1:1:b_m_fnd(mz)
                                        jit = j_idx_divisor(m1,divisor_idx,mz);
                                        ovlp(m1) = junction_info(jit,5)-junction_info(jit,14);
                                    end
                                    [mxv(1), mxi(1)] = max(ovlp);
                                    ovlp(mxi(1)) = 0;
                                    [mxv(2), mxi(2)] = max(ovlp);
                                    if mxv(1) > mxv(2) || mxv(2) == 0
                                        ji = j_idx_divisor(mxi(1),divisor_idx,mz);
                                        ji_tmp = zeros(1,14);
                                        if mz == 1
                                            [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji,:), ji_tmp, b_inter_disp );
                                        else
                                            [b_comp, cm_lst] = f04_compatibility_check( ji_tmp, junction_info(ji,:), b_inter_disp );
                                        end
                                        if sum( b_comp ) > 0
                                            for mx = 1:1:b_m_fnd(mz)
                                                if mx ~= mxi(1)
                                                    jit = j_idx_divisor(mx,divisor_idx,mz);
                                                    junction_info(jit,13) = 0;
                                                end
                                            end
                                            if isempty( cm_lst )
                                            else
                                                seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                    [dividend_idx divisor_idx cm_lst(1) 0 mxv(1) ji];
                                                seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                    [divisor_idx dividend_idx cm_lst(2) 0 mxv(1) ji];
                                            end
                                            b_m_fnd(mz) = 1;
                                        end
                                    else
                                        ji_tmp = zeros(1,14);
                                        if mz == 1
                                            if mxi(1) == mxi(2)
                                                b_m_fnd(mz) = 1;
                                            else
                                                b_m_fnd(mz) = 1; 
                                            end
                                            for my = 1:1:b_m_fnd(mz)
                                                ji = j_idx_divisor(mxi(my),divisor_idx,mz);
                                                [b_comp, cm_lst] = f04_compatibility_check( junction_info(ji,:), ji_tmp, b_inter_disp );
                                                if sum( b_comp ) > 0
                                                    for mx = 1:1:b_m_fnd(mz)
                                                        if mx ~= mxi(my)
                                                            jit = j_idx_divisor(mx,divisor_idx,mz);
                                                            junction_info(jit,13) = 0;
                                                        end
                                                    end
                                                    if isempty( cm_lst )
                                                    else
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(my) ji];
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(my) ji];
                                                    end
                                                else
                                                end
                                            end
                                        else % if mz == 2
                                            if mxi(1) == mxi(2)
                                                b_m_fnd(mz) = 1;
                                            else
                                                b_m_fnd(mz) = 1; 
                                            end
                                            for my = 1:1:b_m_fnd(mz)
                                                ji = j_idx_divisor(mxi(my),divisor_idx,mz);
                                                [b_comp, cm_lst] = f04_compatibility_check( ji_tmp, junction_info(ji,:), b_inter_disp );
                                                if sum( b_comp ) > 0
                                                    for mx = 1:1:b_m_fnd(mz)
                                                        if mx ~= mxi(my)
                                                            jit = j_idx_divisor(mx,divisor_idx,mz);
                                                            junction_info(jit,13) = 0;
                                                        end
                                                    end
                                                    if isempty( cm_lst )
                                                    else
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [dividend_idx divisor_idx cm_lst(1) 0 mxv(my) ji]; 
                                                        seg_connection_mode_cnt = seg_connection_mode_cnt + 1;
                                                        seg_connection_mode_lst( seg_connection_mode_cnt, 1:6) = ...
                                                            [divisor_idx dividend_idx cm_lst(2) 0 mxv(my) ji]; 
                                                    end
                                                end
                                            end
                                        end
                                    end
                                else
                                    b_m_fnd(mz) = 0;
                                end
                            end
                        end
                    end % if b_m_fnd(1) > 1 || b_m_fnd(2) > 1 %n_js > 2
                    if b_m_fnd(1) > 1 || b_m_fnd(2) > 1 %n_js > 2
                        if b_inter_disp > 0
                            str_disp = sprintf('      WARNING: (%d,%d) junctions for (%d, %d)  ', b_m_fnd(1), b_m_fnd(2), dividend_idx, divisor_idx );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            Nchar = 0;
                        end
                    end
                    n_possible_incompatible_junctions = n_possible_incompatible_junctions + max( b_m_fnd(1)-1, 0 ) + max( b_m_fnd(2)-1, 0 );        
                end
                % if mod(n, n_step) == 0
                    % fprintf('.');
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = fprintf( '%s %d/%d/%d (%d/%d)', str_hdr, n, n_seg, round(n_seg_sel), lp, N_Loops ); 
                % end
            end
        end % if seg(dividend_idx).len    
    end % for n = 1:1:n_seg
end % for lp
end
clear j_cnt_divisor;
clear j_idx_divisor;
clear Seg_Nmer_map;
clear Seg_Nmer_seq;
seg_connection_mode_lst = seg_connection_mode_lst(1:seg_connection_mode_cnt, :); 

%% Check group

% set cntg group
% cntg_gid = zeros( cntg_cnt, 1 );
% for k = 1:1:n_seg
%     cntg_gid( seg(k).id ) = group_id(k);
% end
% cntg_gid( cntg_valid_ind < 0 ) = 0;

%% Update cntg_valid_ind
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) < 0
        % cntg_gid(k) = 0;
    else
        if cntg_valid_ind(k) == 0 || cntg_valid_ind(k) == k
            cntg_valid_ind(k) = k;
            % if cntg_gid(k) == 0
            %     fprintf('\n   ERROR1: cntg_valid_ind(%d) == 0 and cntg_gid(%d) == 0 ', k, k );
            % end
        else
            cid_p = k;
            gid_tmp = 0;
            while(1)
                cid = cntg_valid_ind(cid_p);
                if cntg_valid_ind(cid) < 0
                    gid_tmp = -1;
                    break;
                else
                    if cntg_valid_ind(cid) == 0 || cntg_valid_ind(cid) == cid
                        % if cntg_gid(cid) == 0
                        %    fprintf('\n   ERROR2: cntg_valid_ind(%d) == 0 and cntg_gid(%d) == 0 ', cid, cid );
                        % end
                        gid_tmp = cid;
                        break;
                    else
                        cid_p = cid;
                    end
                end
            end
            cid_p = k;
            while(1)
                % cntg_gid(cid_p) = gid_tmp;
                cid = cntg_valid_ind(cid_p);
                cntg_valid_ind(cid_p) = gid_tmp;
                if cntg_valid_ind(cid) <= 0 || cntg_valid_ind(cid) == cid
                    cntg_valid_ind(cid) = gid_tmp;
                    break;
                else
                    cid_p = cid;
                end
            end
        end
    end
end

%% set cid_to_seg_map using cntg_valid_ind
cid_to_seg_map = zeros( cntg_cnt, 1, 'uint32' );
for k = 1:1:n_seg
    cid_to_seg_map( seg(k).id ) = k;
end
N_error = 0;
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) > 0 % && cntg_valid_ind(k) ~= k
        if cid_to_seg_map( k ) == 0
            cid_to_seg_map( k ) = cid_to_seg_map( cntg_valid_ind(k) );
            if cntg_valid_ind( cntg_valid_ind(k) ) ~= cntg_valid_ind(k)
                N_error = N_error + 1;
            end
        else
            if cid_to_seg_map( k ) ~= cid_to_seg_map( cntg_valid_ind(k) )
                N_error = N_error + 1;
            end
        end
    end
end
if N_error > 0
    fprintf(' ERROR: set cid_to_seg_map .. %d ', N_error );
end

%% Set list of possible junctions - Use read_cntg_map, cntg_valid_ind, cid_to_seg_map
Max_n_junctions_per_cntg = 200;
jpc_tmp.bs = Max_n_junctions_per_cntg;
jpc_tmp.nj = 0;
jpc_tmp.ji = zeros( Max_n_junctions_per_cntg, 1, 'uint32' );
jpc = repmat( jpc_tmp, n_seg, 1 );
n_pj = 0;
% if R_mode > 0
Nrt = round( N_reads_total/(R_mode+1) );
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 && read_cntg_map(k,2) > 0
        if read_cntg_map(k,1) ~= read_cntg_map(k,2)
            cntg1 = cntg_valid_ind( read_cntg_map(k,1) );
            cntg2 = cntg_valid_ind( read_cntg_map(k,2) );
            if cntg1 == 0
                cntg1 = read_cntg_map(k,1);
            end
            if cntg2 == 0
                cntg2 = read_cntg_map(k,2);
            end
            if cntg1 > 0 && cntg2 > 0
                if cntg1 ~= cntg2
                    for m = 1:1:2
                        if m == 1
                            cid1 = cid_to_seg_map( cntg1 );
                            cid2 = cid_to_seg_map( cntg2 );
                        else
                            cid2 = cid_to_seg_map( cntg1 );
                            cid1 = cid_to_seg_map( cntg2 );
                        end
                        b_tmp = 1;
                        if jpc( cid1 ).nj > 0
                            jmatch = find( jpc( cid1 ).ji(1:jpc( cid1 ).nj) == cid2, 1 );
                            if ~isempty( jmatch )
                                b_tmp = 0;
                            end
                        end
                        if b_tmp > 0
                            jpc( cid1 ).nj = jpc( cid1 ).nj + 1;
                            if jpc( cid1 ).nj > jpc( cid1 ).bs
                                jpc( cid1 ).ji = [jpc( cid1 ).ji; zeros( Max_n_junctions_per_cntg, 1, 'uint32' ) ];
                                jpc( cid1 ).bs = jpc( cid1 ).bs + Max_n_junctions_per_cntg;
                            end
                            jpc( cid1 ).ji( jpc( cid1 ).nj ) = cid2;
                            if m == 1
                                n_pj = n_pj +1;
                            end
                        end
                    end
                end
            end
        end
    end
end

seg_cnt = 0;
seg_wrong = 0;
for k = 1:1:cntg_cnt
    if cntg_valid_ind(k) == k 
        seg_cnt = seg_cnt +1;
    end
    if cntg_valid_ind(k) == 0 
        seg_wrong = seg_wrong + 1;
    end
end

% end
%% Check group

b_group_valid = zeros(n_seg,1);
for n = 1:1:n_seg
    m = group_id(n);
    b_group_valid(m) = 1;
end

if N_groups == sum( b_group_valid )
    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id(k);
        Group_Size(m) = Group_Size(m) + 1;
    end
    str_disp = sprintf('done, N_groups/N_junctions: %d(%d)/%d ', ...
        N_groups, max(Group_Size), N_junctions );
    fprintf(repmat('\b', 1, Nchar));
    fprintf( '%s %s', str_hdr, str_disp ); 
    fprintf( fp_log, ' %s', str_disp ); 
else
    N_groups_tmp = 0;
    for n = 1:1:N_groups
        if b_group_valid(n) == 0 
        else
            N_groups_tmp = N_groups_tmp + 1;
            if N_groups_tmp ~= n
                jidxs = find( junction_info(1:N_junctions,1) == n );
                if ~isempty(jidxs)
                    junction_info(jidxs,1) = N_groups_tmp;
                end
                sidxs = find( group_id(1:n_seg) == n );
                if ~isempty(sidxs)
                    group_id(sidxs) = N_groups_tmp;
                end
                b_group_valid(N_groups_tmp) = 1;
                b_group_valid(n) = 0;
            end
        end
    end
    N_groups = N_groups_tmp;
    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id(k);
        Group_Size(m) = Group_Size(m) + 1;
    end
    str_disp = sprintf('done, N_groups/N_junctions: %d(%d)/%d ', ...
        N_groups, max(Group_Size), N_junctions );
    fprintf(repmat('\b', 1, Nchar));
    fprintf( '%s %s', str_hdr, str_disp ); 
    fprintf( fp_log, ' %s', str_disp ); 
end
fprintf(' Npj: %d, (%d, %d, %d) ', n_pj, n_seg, seg_cnt, seg_wrong );

%% Generate read_group_map

Nrt = round( N_reads_total/(R_mode+1) );
Nrt_new = 0;
Nrt_mapped = 0;
Nrt_mapped1 = 0;
Nrt_wrong = 0;
Nrt_gs = 0;

% fpx = fopen('test.txt', 'wt' );
kmod = round(Nrt/6);
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 || read_cntg_map(k,2) > 0
        Nrt_mapped = Nrt_mapped + 1;
        if read_cntg_map(k,1) <= 0 || read_cntg_map(k,2) <= 0
            Nrt_mapped1 = Nrt_mapped1 + 1;
        end
        if read_cntg_map(k,1) > 0 && read_cntg_map(k,2) > 0
            if cid_to_seg_map( read_cntg_map(k,1) ) > 0 && cid_to_seg_map( read_cntg_map(k,2) ) > 0
                if group_id( cid_to_seg_map( read_cntg_map(k,1) ) ) ~= group_id( cid_to_seg_map( read_cntg_map(k,2) ) )
                    Nrt_wrong = Nrt_wrong + 1;
%                     fprintf( fpx, '%8d (%f, %d)\t%8d (%f, %d)\n', cid_to_seg_map( read_cntg_map(k,1) ), ...
%                         seg(cid_to_seg_map( read_cntg_map(k,1) )).ave_cvg_dep, seg(cid_to_seg_map( read_cntg_map(k,1) )).len, ...
%                         cid_to_seg_map( read_cntg_map(k,2) ), ...
%                         seg(cid_to_seg_map( read_cntg_map(k,2) )).ave_cvg_dep, seg(cid_to_seg_map( read_cntg_map(k,2) )).len );
                else
                    Nrt_new = Nrt_new + 1;
                    if Group_Size( group_id( cid_to_seg_map( read_cntg_map(k,1) ) ) ) > 1
                        Nrt_gs = Nrt_gs + 1;
                    end
                end
            end
        end
    end
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
% fclose(fpx);

% fpx = fopen('test_j.txt', 'wt' );
% for k = 1:1:n_seg
%     fprintf(fpx, '%d - ', k );
%     for m = 1:1:jpc(k).nj
%         fprintf(fpx, '%d\t', jpc(k).ji(m) );
%     end
%     fprintf(fpx, '\n' );
% end
% fclose(fpx);
fprintf(' done. ');

%% Filter incompatible junctions

% n_jcnt = 0;
% for n = 1:1:N_junctions
%     if junction_info(n,13) > 0
%         n_jcnt = n_jcnt + 1;
%     end
% end
n_jcnt = sum( junction_info(1:N_junctions, 13) > 0 );

str_disp = sprintf('\n   Invalidating .. %d(%d)/%d ', N_groups, max(Group_Size), n_jcnt );

fprintf( '%s', str_disp ); 
fprintf( fp_log, ' %s', str_disp ); 

[ji_sorted, jidx_sorted] = sort( junction_info(1:N_junctions,3), 'ascend' ); 
ptr_divisors = zeros( n_seg, 3, 'int32' );
cur_ptr = 1;
ptr_divisors( cur_ptr, 1 ) = jidx_sorted(1);
ptr_divisors( cur_ptr, 2 ) = 1;
for n = 2:1:N_junctions
    if junction_info(jidx_sorted(n),3) ~= cur_ptr
        ptr_divisors( cur_ptr, 3 ) = n-1;
        cur_ptr = cur_ptr + 1;
        ptr_divisors( cur_ptr, 1 ) = jidx_sorted(n);
        ptr_divisors( cur_ptr, 2 ) = n;
    end
    if n == N_junctions
        ptr_divisors( cur_ptr, 3 ) = n;
    end
end
N_divisors = cur_ptr;
for n = 1:1:N_divisors
    if ptr_divisors( cur_ptr, 2 ) < ptr_divisors( cur_ptr, 3 )
        sp = ptr_divisors( cur_ptr, 2 );
        ep = ptr_divisors( cur_ptr, 3 );
        nj = ep-sp+1;
        msk = zeros(nj,2, class(junction_info) );
        for k = 1:1:nj
            jm = junction_info(jidx_sorted(sp+k-1),12);
            if jm == 1 || jm == 3 % head
                msk(k,1) = 1;
            else % tail
                msk(k,2) = 1;
            end
        end
        Ovlp_len = junction_info(jidx_sorted(sp:ep),5) - junction_info(jidx_sorted(sp:ep),14);
        for m = 1:1:2
            if sum( msk(:,m) ) > 1
                [mxv,mxi] = max( Ovlp_len.*msk(:,m) );
                for k = 1:1:nj
                    jm = junction_info(jidx_sorted(sp+k-1),12);
                    if jm == m || jm == (m+2)
                        if k ~= mxi
                            junction_info(jidx_sorted(sp+k-1),13) = 0;
                        end
                    end
                end
            end
        end
    end
end

% n_jcnt = 0;
% for n = 1:1:N_junctions
%     if junction_info(n,13) > 0
%         n_jcnt = n_jcnt + 1;
%     end
% end
n_jcnt = sum( junction_info(1:N_junctions, 13) > 0 );
str_disp = sprintf('-> %d(%d)/%d ', N_groups, max(Group_Size), n_jcnt );
fprintf( '%s', str_disp ); 
fprintf( fp_log, '%s', str_disp ); 

fprintf('-');    
while(1)
    % regroup
    group_id = zeros(n_seg, 1);
    ndd = zeros(n_seg, 2);
    junction_info(:,1) = 0;
    b_group_valid = ones( n_seg, 1, 'int8' );
    for n = 1:1:N_junctions
        if junction_info(n,13) > 0
            ndd( junction_info(n,2), 1 ) = ndd( junction_info(n,2), 1 ) + 1;
            ndd( junction_info(n,3), 2 ) = ndd( junction_info(n,3), 2 ) + 1;
            if group_id( junction_info(n,2) ) == 0
                if group_id( junction_info(n,3) ) == 0
                    N_groups_tmp = find( b_group_valid, 1 );
                    b_group_valid( N_groups_tmp ) = 0;
                    group_id( junction_info(n,2) ) = N_groups_tmp;
                    group_id( junction_info(n,3) ) = N_groups_tmp;
                    junction_info(n,1) = N_groups_tmp;
                else
                    group_id( junction_info(n,2) ) = group_id( junction_info(n,3) );
                    junction_info(n,1) = group_id( junction_info(n,3) );
                end
            else
                if group_id( junction_info(n,3) ) == 0
                    junction_info(n,1) = group_id( junction_info(n,2) );
                    group_id( junction_info(n,3) ) = junction_info(n,1);
                else
                    % ERROR
                    if group_id( junction_info(n,3) ) == group_id( junction_info(n,2) )
                        junction_info(n,1) = group_id( junction_info(n,2) );
                    else
                        junction_info(n,1) = group_id( junction_info(n,2) );
                        gid_tmp = group_id( junction_info(n,3) );
                        for m = 1:1:n
                            if junction_info(m,1) == gid_tmp
                                junction_info(m,1) = junction_info(n,1);
                            end
                        end
                        for m = 1:1:n_seg
                            if group_id(m) == gid_tmp
                                group_id(m) = junction_info(n,1);
                            end
                        end
                    end
                end
            end
        end
    end

    for n = 1:1:n_seg
        if group_id(n) == 0
            N_groups_tmp = find( b_group_valid, 1 );
            b_group_valid( N_groups_tmp ) = 0;
            group_id(n) = N_groups_tmp;
        end
    end
    N_groups = sum( 1 - b_group_valid );
    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id(k);
        Group_Size(m) = Group_Size(m) + 1;
    end
    nkp = 0;
    if R_mode >= 0
        if max( Group_Size ) <= cfg.max_n_cntg_per_group
            break;
        else
            gidx = find( Group_Size > cfg.max_n_cntg_per_group );
            % ovlp_chk = ones( length(gidx), 1 );
            b_tmp = 0;
            for k = 1:1:length(gidx)
                jidx = find( junction_info(:,1) == gidx(k) );
                % unsel = (ndd( junction_info(jidx,3), 1) == 0).*(ndd( junction_info(jidx,3), 2 ) == 1);
                % if sum( unsel ) < length(unsel)
                    % min_ovlp = min( junction_info(jidx,5) - junction_info(jidx,14) ) + int32(10000.*int32(unsel)) );
                    min_ovlp = min( junction_info(jidx,5) - junction_info(jidx,14) ); 
                    % max_dist = max( junction_info(jidx,14) );
                    if min_ovlp < cfg.safe_overlap_threshold || Group_Size(gidx(k)) > round(cfg.max_n_cntg_per_group*1.33)
                        for m = 1:1:length(jidx)
                            if junction_info(jidx(m),5) - junction_info(jidx(m),14) <= min_ovlp
                            % if junction_info(jidx(m),14) >= max_dist
                                % if ndd( junction_info(jidx(m),3), 1 ) == 0 && ndd( junction_info(jidx(m),3), 2 ) == 1
                                % else
                                    junction_info(jidx(m),13) = 0;
                                    junction_info(jidx(m),1) = 0;
                                    b_tmp = 1;
                                % end
                            end
                        end
                    end
                % else
                % end
            end
            if b_tmp == 0
                break;
            else
                fprintf('>');
            end
        end
    else
        if max( Group_Size ) <= cfg.max_n_cntg_per_group
            break;
        else
            gidx = find( Group_Size > cfg.max_n_cntg_per_group );
            % ovlp_chk = ones( length(gidx), 1 );
            b_tmp = 0;
            for k = 1:1:length(gidx)
                jidx = find( junction_info(:,1) == gidx(k) );
                % unsel = (ndd( junction_info(jidx,3), 1) == 0).*(ndd( junction_info(jidx,3), 2 ) == 1);
                % if sum( unsel ) < length(unsel)
                    % min_ovlp = min( junction_info(jidx,5) - junction_info(jidx,14) ) + int32(10000.*int32(unsel)) );
                    min_ovlp = min( junction_info(jidx,5) - junction_info(jidx,14) ); 
                    % max_dist = max( junction_info(jidx,14) );
                    if min_ovlp < cfg.safe_overlap_threshold || Group_Size(gidx(k)) > round(cfg.max_n_cntg_per_group*1.33)
                        for m = 1:1:length(jidx)
                            if junction_info(jidx(m),5) - junction_info(jidx(m),14) <= min_ovlp
                            % if junction_info(jidx(m),14) >= max_dist
                                nj = jpc( junction_info(jidx(m),2) ).nj;
                                if nj == 0
                                    junction_info(jidx(m),13) = 0;
                                    junction_info(jidx(m),1) = 0;
                                    b_tmp = 1;
                                else
                                    if isempty( find( jpc( junction_info(jidx(m),2) ).ji(1:nj) ==  junction_info(jidx(m),3), 1 ) )
                                    % else
                                        junction_info(jidx(m),13) = 0;
                                        junction_info(jidx(m),1) = 0;
                                        b_tmp = 1;
                                    else
                                        nkp = nkp + 1;
                                    end
                                end
                            end
                        end
                    end
                % else
                % end
            end
            if b_tmp == 0
                break;
            else
                fprintf('>');
            end
        end
    end
end

% n_jcnt = 0;
% for n = 1:1:N_junctions
%     if junction_info(n,13) > 0
%         n_jcnt = n_jcnt + 1;
%     end
% end
n_jcnt = sum( junction_info(1:N_junctions, 13) > 0 );

% str_disp = sprintf('> %d(%d)/%d (MO: %d/%d) %d ', N_groups, max(Group_Size), n_jcnt, min_ovlp, max_ovlp, nkp );
str_disp = sprintf('> %d(%d)/%d (MO: %d/%d) ', N_groups, max(Group_Size), n_jcnt, min_ovlp, max_ovlp );
fprintf( '%s', str_disp ); 
fprintf( fp_log, '%s', str_disp ); 

scm_cnt = 0;
for n = 1:1:seg_connection_mode_cnt
    ji = seg_connection_mode_lst(n,6);
    if junction_info(ji,13) ~= 0
        scm_cnt = scm_cnt + 1;
        if scm_cnt ~= n
            seg_connection_mode_lst(scm_cnt,:) = seg_connection_mode_lst(n,:);
        end
    end
end
seg_connection_mode_cnt = scm_cnt;
seg_connection_mode_lst = seg_connection_mode_lst(1:seg_connection_mode_cnt, :); 

%% Strand correction using junction info.
fprintf( '\n   Checking compatibility' );
str_hdr = sprintf(' .....' );
Nchar = fprintf( '%s', str_hdr );
fprintf( fp_log, '\n   Checking compatibility%s', str_hdr );
% input
%seg_connection_mode = seg_connection_mode + seg_connection_mode';
%seg_connection_depth = seg_connection_depth + seg_connection_depth';
% Outputput
Strand_corrected = int8( zeros(n_seg, 1) );
% Group_index = zeros(n_seg, 1);

% Temporary
seg_connection_mode_lst(1:seg_connection_mode_cnt,4) = group_id( seg_connection_mode_lst(1:seg_connection_mode_cnt,1) );
grp_chk = group_id( seg_connection_mode_lst(1:seg_connection_mode_cnt,2) );
if sum( seg_connection_mode_lst(1:seg_connection_mode_cnt,4) ~= grp_chk ) > 0
    fprintf('\n   Group check error ');
end
clear grp_chk;
sum_chk = sum( seg_connection_mode_lst(1:seg_connection_mode_cnt,4) == 0 );
if sum_chk > 0
    fprintf('\n   Sum check error ');
    d_cnt = 0;
    for k = 1:1:seg_connection_mode_cnt
        if seg_connection_mode_lst(k,4) > 0
            d_cnt = d_cnt + 1;
            seg_connection_mode_lst(d_cnt,:) = seg_connection_mode_lst(k,:);
        end
    end
    seg_connection_mode_cnt = d_cnt;
end
clear sum_chk;
[axx, sorted_gidx] = sort( seg_connection_mode_lst(1:seg_connection_mode_cnt,4), 'descend' );
clear axx;

idx_ptr = 1;
gidx_map = zeros( n_seg, 1, 'uint32' );
% n_step = round(N_groups/8);
Num_trial = 1;
for n = 1:1:N_groups
    gid_tmp = seg_connection_mode_lst( sorted_gidx( idx_ptr ), 4 );
    n_cnt = 0;
    while(1)
        if idx_ptr+n_cnt > seg_connection_mode_cnt
            break;
        else
            if seg_connection_mode_lst( sorted_gidx( idx_ptr+n_cnt ), 4 ) == gid_tmp 
                n_cnt = n_cnt + 1;
            else
                break;
            end
        end
    end
    if n_cnt == 0
        break;
    end
    if Group_Size(gid_tmp) == 1
        fprintf('\n    ERROR: Group_Size(gid_tmp) == 1' );
    else
        seg_idx = find( group_id == gid_tmp );
        if length( seg_idx ) ~= Group_Size(gid_tmp)
            fprintf('\n    ERROR:  length( seg_idx ) ~= Group_Size(gid_tmp)' );
        end
        seg_connection_mode_selected = zeros(Group_Size(gid_tmp), Group_Size(gid_tmp), 'int8' ); 
        % ncnt = 0;
        for m = 1:1:n_cnt %Group_Size(gid_tmp)
            sci = seg_connection_mode_lst( sorted_gidx(idx_ptr+m-1), : );
            cidx = find( seg_idx == sci(1) );
            if gidx_map( sci(1) ) == 0
                gidx_map( sci(1) ) = cidx;
            else
                if gidx_map( sci(1) ) ~= cidx
                    fprintf('\n    ERROR: gidx_map( sci(1) ) ~= cidx');
                end
            end
        end
        if sum( gidx_map( seg_idx ) == 0 ) > 0
            fprintf('\n    WARNING: sum( gidx_map( seg_idx ) == 0 ) > 0 ');
            b_seg_mapped = zeros( Group_Size(gid_tmp), 1 );
            for m = 1:1:Group_Size(gid_tmp)
                if gidx_map( seg_idx(m) ) > 0 && gidx_map( seg_idx(m) ) <= Group_Size(gid_tmp)
                    b_seg_mapped(m) = 1;
                end
            end
            for m = 1:1:Group_Size(gid_tmp)
                if gidx_map( seg_idx(m) ) == 0
                    cidx = find( b_seg_mapped(m) == 0, 1 );
                    if isempty(cidx)
                        fprintf(' --> Could NOT correct ERROR ');
                        break;
                    else
                        gidx_map( seg_idx(m) ) = cidx;
                        b_seg_mapped(m) = 1;
                    end
                    if sum( b_seg_mapped ) == Group_Size(gid_tmp)
                        break;
                    end
                end
            end
        end
        for m = 1:1:n_cnt 
            sci = seg_connection_mode_lst( sorted_gidx(idx_ptr+m-1), : );
            seg_connection_mode_selected( gidx_map(sci(1)), gidx_map(sci(2)) ) = sci(3);
        end
        d_flag = 0;
        ovlp = seg_connection_mode_lst( sorted_gidx(idx_ptr:idx_ptr+n_cnt-1), 5 );
        seg_connection_mode_selected_save = seg_connection_mode_selected;
        rem_idx = 0;
        N_loops = min( Num_trial, n_cnt+1 );
        for kk = 1:N_loops
            seg_connection_mode_selected = seg_connection_mode_selected_save;
            if rem_idx > 0 && kk < N_loops
                sci2 = seg_connection_mode_lst( sorted_gidx(rem_idx), : );
                seg_connection_mode_selected( gidx_map(sci2(1)), gidx_map(sci2(2)) ) = 0;
                seg_connection_mode_selected( gidx_map(sci2(2)), gidx_map(sci2(1)) ) = 0;
            else
                rem_idx = 0;
            end
            c_tmp = 0;
            Que_seg_idx = zeros(n_seg, 1,'uint32');
            Q_len = 0;
            Strand_corrected(seg_idx(1:Group_Size(gid_tmp))) = 0;
            while(1)
                if Q_len == 0
                    % get a seg index
                    b_tmp = 0;
                    for k = 1:1:Group_Size(gid_tmp)
                        m = seg_idx(k);
                        if Strand_corrected(m) == 0
                            b_tmp = 1;
                            break;
                        end
                    end
                    if b_tmp == 0
                        break;
                    else
                        ref_idx = m;
                        Strand_corrected(ref_idx) = 1;
                        for k = 1:1:Group_Size(gid_tmp)
                            m = seg_idx(k);
                            if seg_connection_mode_selected(gidx_map(ref_idx), gidx_map(m)) ~= 0
                                [axx, strand_tmp] = ...
                                    f04_strand_correction( Strand_corrected(ref_idx), seg_connection_mode_selected(gidx_map(ref_idx), gidx_map(m)) );
                                if Strand_corrected(m) == 0
                                    Q_len = Q_len + 1;
                                    Que_seg_idx( Q_len ) = m;
                                    Strand_corrected(m) = strand_tmp;
                                else
                                    if Strand_corrected(m) ~= strand_tmp 
                                        if b_inter_disp > 0 && d_flag == 0
                                            d_flag = 1;
                                            str_disp = sprintf( '      ERROR: Direction ambiguous in Group %d, Cnt: %d', group_id(ref_idx), kk );
                                            fprintf('\n%s', str_disp );
                                        end
                                        c_tmp = 1;
                                        if kk < N_loops 
                                            break;
                                        end
                                    else
                                    end
                                end
                            end
                        end
                        if c_tmp == 1 && kk < N_loops 
                            break;
                        end
                    end
                else % if Queue is not empty
                    ref_idx = Que_seg_idx(1);
                    Que_seg_idx(1:Q_len-1) = Que_seg_idx(2:Q_len);
                    Q_len = Q_len - 1;

                    for k = 1:1:Group_Size(gid_tmp)
                        m = seg_idx(k);
                        if seg_connection_mode_selected(gidx_map(ref_idx), gidx_map(m)) ~= 0
                            [axx, strand_tmp] = ...
                                f04_strand_correction( Strand_corrected(ref_idx), seg_connection_mode_selected(gidx_map(ref_idx), gidx_map(m)) );
                            if Strand_corrected(m) == 0
                                Q_len = Q_len + 1;
                                Que_seg_idx( Q_len ) = m;
                                Strand_corrected(m) = strand_tmp;
                            else
                                if Strand_corrected(m) ~= strand_tmp 
                                    if b_inter_disp > 0 && d_flag == 0
                                        d_flag = 1;
                                        str_disp = sprintf( '      ERROR: Direction ambiguous in Group %d, Cnt: %d', group_id(ref_idx), kk );
                                        fprintf('\n%s', str_disp );
                                    end
                                    c_tmp = 1;
                                    if kk < N_loops 
                                        break;
                                    end
                                else
                                end
                            end
                        end
                    end
                    if c_tmp > 0 && kk < N_loops 
                        break;
                    end
                end
            end
            if c_tmp > 0
                [axx, mni] = min( ovlp );
                ovlp(mni) = 100000;
                rem_idx = mni;
                % str_disp = sprintf( '      WARNING: Incompatible connections detected in Group %d (%d/%d) ', group_id(ref_idx), kk, n_cnt );
                % fprintf('\n%s', str_disp );
            else
                % str_disp = sprintf( '      WARNING: Incompatible connections corrected in Group %d (%d/%d) ', group_id(ref_idx), kk, n_cnt );
                % fprintf('\n%s', str_disp );
                if rem_idx > 0
                    j_idx = seg_connection_mode_lst( sorted_gidx(rem_idx), 6 );
                    junction_info(j_idx,13) = 0;
                end
                break;
            end
        end
        if c_tmp > 0 && Num_trial > 1
            str_disp = sprintf( '      WARNING: Couldn''t resolve Incompatible junctions in Group %d  ', group_id(ref_idx) );
            fprintf('\n%s', str_disp );
            fprintf(fp_log, '\n%s', str_disp );
            Nchar = 0;
        end
    end
    idx_ptr = idx_ptr + n_cnt; 
    if idx_ptr > seg_connection_mode_cnt
        fprintf(repmat('\b', 1, Nchar));
        fprintf( '%s %d/%d ', str_hdr, N_groups, N_groups ); 
        break;
    else
        fprintf(repmat('\b', 1, Nchar));
        Nchar = fprintf( '%s %d/%d ', str_hdr, n, N_groups ); 
    end
    % if mod(n, n_step) == 0
        % fprintf('.');
    % end
end
for k = 1:1:n_seg
    if Strand_corrected(k) == 0
        Strand_corrected(k) = 1;
    end
end
clear Que_seg_idx;

% Correct strand
for k = 1:1:n_seg
    seg_idx = k;
    if Strand_corrected(seg_idx) < 0
        seg(seg_idx).cvg_dep(1:4,1:seg(seg_idx).len) = seg(seg_idx).cvg_dep(4:-1:1,seg(seg_idx).len:-1:1);
        seg(seg_idx).seq(1:seg(seg_idx).len) = f03_seq_est( seg(seg_idx).cvg_dep(:,1:seg(seg_idx).len) );
        seg(seg_idx).ave_cvg_dep = mean( sum( seg(seg_idx).cvg_dep(:,1:seg(seg_idx).len) ) );
    end
end

%% Update junction info. according to Strand correction
n_step = round(N_junctions/8);
D_flag = zeros( N_groups, 1, 'int8' );
N_junctions_eff = 0;
for k = 1:1:N_junctions
    if junction_info(k,13) > 0
        N_junctions_eff = N_junctions_eff +1;
        grp = junction_info(k,1);
        dividend_idx = junction_info(k,2);
        divisor_idx = junction_info(k,3);
        o_len = junction_info(k,5);
        j_mod = junction_info(k,12);

        % junction_mode = 1 -> divisor head, split
        % junction_mode = 2 -> divisor tail, merge 
        % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
        % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
        if Strand_corrected(dividend_idx) > 0
            % No action
            switch( j_mod )
                case 1,
                    if Strand_corrected(divisor_idx) > 0
                        % No action
                    else
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    end
                case 2,
                    if Strand_corrected(divisor_idx) > 0
                        % No action
                    else
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    end
                case 3,
                    if Strand_corrected(divisor_idx) > 0
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    else
                        junction_info(k,9) = seg(divisor_idx).len - o_len + 1;
                        junction_info(k,10) = seg(divisor_idx).len;
                    end
                otherwise,
                    if Strand_corrected(divisor_idx) > 0
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    else
                        junction_info(k,9) = 1;
                        junction_info(k,10) = o_len;
                    end
            end
        else
            pos = junction_info(k,6);
            junction_info(k,6) = seg(dividend_idx).len - (pos+o_len-1) + 1;
            junction_info(k,7) = seg(dividend_idx).len - pos + 1;

            % junction_mode = 1 -> divisor head, split
            % junction_mode = 2 -> divisor tail, merge 
            % junction_mode = 3 -> divisor head, merge -> reverse divisor strand
            % junction_mode = 4 -> divisor tail, split -> reverse divisor strand 
            switch( j_mod )
                case 1,
                    if Strand_corrected(divisor_idx) > 0
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                             % Nchar = 0;
                       end
                    else
                        junction_info(k,9) = seg(divisor_idx).len - o_len + 1;
                        junction_info(k,10) = seg(divisor_idx).len;
                    end
                case 2,
                    if Strand_corrected(divisor_idx) > 0
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                             % Nchar = 0;
                       end
                    else
                        junction_info(k,9) = 1;
                        junction_info(k,10) = o_len;
                    end
                case 3,
                    if Strand_corrected(divisor_idx) > 0
                        if junction_info(k,9) == 1 && junction_info(k,10) == o_len
                            % OK
                        else
                            if b_inter_disp > 0 && D_flag(grp) == 0
                                D_flag(grp) = 1;
                                str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                                fprintf('\n%s', str_disp );
                                fprintf(fp_log, '\n%s', str_disp );
                                % Nchar = 0;
                            end
                        end
                    else
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    end
                otherwise,
                    if Strand_corrected(divisor_idx) > 0
                        if junction_info(k,9) == (seg(divisor_idx).len-o_len+1) && junction_info(k,10) == seg(divisor_idx).len
                            % OK
                        else
                            if b_inter_disp > 0 && D_flag(grp) == 0
                                D_flag(grp) = 1;
                                str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                                fprintf('\n%s', str_disp );
                                fprintf(fp_log, '\n%s', str_disp );
                                % Nchar = 0;
                            end
                        end
                    else
                        if b_inter_disp > 0 && D_flag(grp) == 0
                            D_flag(grp) = 1;
                            str_disp = sprintf( '      WARNING: Incompatible Junctions in Group %d  ', grp );
                            fprintf('\n%s', str_disp );
                            fprintf(fp_log, '\n%s', str_disp );
                            % Nchar = 0;
                        end
                    end
            end
        end
        if mod(k, n_step) == 0
            fprintf('.');
            fprintf(fp_log, '.');
        end
    end
end
clear D_flag;
%fprintf('%s', 'Strand correction completed ...' );
Junction_info = junction_info(1:N_junctions,:);
 
%% Checking integrity
% str_disp = sprintf(' Integrity check ... ' );
% fprintf( '%s', str_disp );
% fprintf( fp_log, '%s', str_disp );

N_wrong = 0;
Max_diff = 0;
N_junctions_eff = 0;
% Dth = cfg.norm_dist_threshold*cfg.nominal_read_length;
for k = 1:1:N_junctions
    % get junction info
    if Junction_info(k,13) > 0
        N_junctions_eff = N_junctions_eff + 1;
        ref_idx = ( Junction_info(k,2) );
        tgt_idx = ( Junction_info(k,3) );
        jlen = Junction_info(k,5);
        ref_pos = Junction_info(k,6);
        tgt_pos = Junction_info(k,9);

        sr = seg(ref_idx).seq(ref_pos:1:ref_pos+jlen-1);
        st = seg(tgt_idx).seq(tgt_pos:1:tgt_pos+jlen-1);
        msk = ( abs(sign( sr - st )) );
        smsk = sum(msk);
        Dth = cfg.norm_dist_threshold*jlen;
        % init junction_mat
        if smsk > Dth
            Junction_info(k,13) = 0;
            N_wrong = N_wrong + 1;
            if smsk > Max_diff
                Max_diff = smsk;
            end
        end
    end
end
if cfg.max_n_cntg_per_group > 0
    [axx, group_idx_sorted] = sort(Group_Size, 'descend');
    N_seg_max_per_group = round(cfg.max_n_cntg_per_group*1.5);
    for n = 1:1:N_groups
        group_idx = group_idx_sorted(n); 
        % get seg indices for this group
        if Group_Size(group_idx) > N_seg_max_per_group
            seg_idxs = zeros(Group_Size(group_idx),1, 'uint32');
            n_seg_in_group = 0;
            for k = 1:1:n_seg
                if group_id(k) == group_idx
                    n_seg_in_group = n_seg_in_group + 1;
                    seg_idxs(n_seg_in_group) = k;
                end
            end
            seg_cvg = [seg(seg_idxs).ave_cvg_dep].*[seg(seg_idxs).len];
            [axx,srt_idxs] = sort(seg_cvg, 'descend');
            srt_idxs_sel = srt_idxs(N_seg_max_per_group+1:end);
            jidxs = find( Junction_info(1:N_junctions,1) == group_idx );
            for k = 1:1:length(jidxs)
                sm1 = sum( seg_idxs(srt_idxs_sel) == Junction_info(jidxs(k),2) );
                sm2 = sum( seg_idxs(srt_idxs_sel) == Junction_info(jidxs(k),3) );
                if sm1 > 0 || sm2 > 0
                    Junction_info(jidxs(k),13) = 0;
                end
            end
            for k = 1:1:length(srt_idxs_sel)
                group_id( seg_idxs(srt_idxs_sel(k)) ) = N_groups + k;
            end
            N_groups = N_groups + length(srt_idxs_sel);
            str_disp = sprintf( '      WARNING: Too many contigs in Group %d: %d -> %d ', group_idx, Group_Size(group_idx), Group_Size(group_idx)-length(srt_idxs_sel) );
            fprintf('\n%s', str_disp );
            fprintf(fp_log, '\n%s', str_disp );
        end
    end
    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id(k);
        Group_Size(m) = Group_Size(m) + 1;
    end
end

Nchar = 0;
fprintf(repmat('\b', 1, Nchar));
str_disp = sprintf(' done (%d grps, %d ejs) ', N_groups, N_junctions_eff );
fprintf( '%s', str_disp );
fprintf( fp_log, '%s', str_disp );

%% Generate read_group_map

Nrt = round( N_reads_total/(R_mode+1) );
Nrt_new = 0;
Nrt_mapped = 0;
Nrt_mapped1 = 0;
Nrt_wrong = 0;
Nrt_gs = 0;

% fpx = fopen('test2.txt', 'wt' );

kmod = round(Nrt/6);
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 || read_cntg_map(k,2) > 0
        Nrt_mapped = Nrt_mapped + 1;
        if read_cntg_map(k,1) <= 0 || read_cntg_map(k,2) <= 0
            Nrt_mapped1 = Nrt_mapped1 + 1;
        end
        if read_cntg_map(k,1) > 0 && read_cntg_map(k,2) > 0
            if cid_to_seg_map( read_cntg_map(k,1) ) > 0 && cid_to_seg_map( read_cntg_map(k,2) ) > 0
                if group_id( cid_to_seg_map( read_cntg_map(k,1) ) ) ~= group_id( cid_to_seg_map( read_cntg_map(k,2) ) )
                    Nrt_wrong = Nrt_wrong + 1;
%                     fprintf( fpx, '%8d (%f, %d)\t%8d (%f, %d)\n', cid_to_seg_map( read_cntg_map(k,1) ), ...
%                         seg(cid_to_seg_map( read_cntg_map(k,1) )).ave_cvg_dep, seg(cid_to_seg_map( read_cntg_map(k,1) )).len, ...
%                         cid_to_seg_map( read_cntg_map(k,2) ), ...
%                         seg(cid_to_seg_map( read_cntg_map(k,2) )).ave_cvg_dep, seg(cid_to_seg_map( read_cntg_map(k,2) )).len );
                else
                    Nrt_new = Nrt_new + 1;
                    if Group_Size( group_id( cid_to_seg_map( read_cntg_map(k,1) ) ) ) > 1
                        Nrt_gs = Nrt_gs + 1;
                    end
                end
            end
        end
    end
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
% fclose(fpx);
Nrt_new = Nrt_mapped;

fname_out2 = sprintf('%s', out_file_prefix );
fname_txt = sprintf('%s.rgmap', fname_out2 );
fp_t = fopen( fname_txt, 'wt' );
fprintf(fp_t, '%d\n', Nrt_new);

kmod = round(Nrt/6);
for k = 1:1:Nrt
    if read_cntg_map(k,1) > 0 || read_cntg_map(k,2) > 0
        fprintf(fp_t, '%d\t', k );
        if read_cntg_map(k,1) > 0
            if cid_to_seg_map( read_cntg_map(k,1) ) > 0
                fprintf(fp_t, '%d\t', group_id( cid_to_seg_map( read_cntg_map(k,1) ) ) );
            else
                fprintf(fp_t, '%d\t', cid_to_seg_map( read_cntg_map(k,1) ) );
            end
        else
            fprintf(fp_t, '%d\t', read_cntg_map(k,1) );
        end
        
        if read_cntg_map(k,2) > 0
            if cid_to_seg_map( read_cntg_map(k,2) ) > 0
                fprintf(fp_t, '%d\n', group_id( cid_to_seg_map( read_cntg_map(k,2) ) ) );
            else
                fprintf(fp_t, '%d\n', cid_to_seg_map( read_cntg_map(k,2) ) );
            end
        else
            fprintf(fp_t, '%d\n', read_cntg_map(k,2) );
        end
        if mod(k,kmod) == 0
            fprintf('.');
        end
    end
end
fclose(fp_t);

str_disp = sprintf(' (%d, %d, %d, %d, %d) ', Nrt, Nrt_mapped, Nrt_mapped1, Nrt_gs, Nrt_wrong );
fprintf( '%s', str_disp );
fprintf( fp_log, '%s', str_disp );

%% Graph construction

Seg_mem_size = n_seg;
seg_n = repmat( sg, Seg_mem_size, 1 );
group_id_new = zeros(Seg_mem_size,1);

fprintf('\n   Constructing graph' );
str_hdr = sprintf(' ..... ' );
Nchar = fprintf( '%s', str_hdr );
fprintf( fp_log, '\n   Constructing graph%s', str_hdr );

[axx, group_idx_sorted] = sort(Group_Size, 'descend');
Group_Size_new = Group_Size;
n_seg_found = 0;
n_seg_added = 0;
seg_connection_mat_lst = zeros( Seg_mem_size, 3, 'uint32' ); 
seg_connection_mat_cnt = 0;
group_cnt_pos = zeros( N_groups, 3, 'uint32' );
sid2cid = zeros(n_seg,1);

% n_step = round(N_groups/8);
for n = 1:1:N_groups
    group_idx = group_idx_sorted(n); % Must not be sorted
    % get seg indices for this group
    if Group_Size(group_idx) > 1
        
    seg_idxs = find( group_id(1:n_seg) == group_idx ); 
    n_seg_in_group = length(seg_idxs);
    for k = 1:1:n_seg_in_group
        sid2cid(seg_idxs(k)) = k;
    end
    max_seg_length = max( [seg(seg_idxs).len] );
    
    % initialize junction_mat
    junction_matrix_tmp.mat = uint32(0);
    junction_matrix = repmat( junction_matrix_tmp, n_seg_in_group, 1 );
    mem_req = 0;
    for k = 1:1:n_seg_in_group
        mem_req = mem_req + n_seg_in_group*seg(seg_idxs(k)).len*4;
    end    
    if mem_req > 10*1000000000
        str = sprintf('   WARNING: %d GBytes memory required ', ceil(mem_req/1000000000) );
        fprintf('\n   %s', str);
        fprintf(fp_log, '\n   %s', str);
    end
    for k = 1:1:n_seg_in_group
        junction_matrix(k).mat = zeros(n_seg_in_group, seg(seg_idxs(k)).len,'uint32');
    end    
    jsel = find( Junction_info(1:N_junctions,1) == group_idx );
    for k1 = 1:1:length(jsel)
        k = jsel(k1);
        if Junction_info(k,1) == group_idx && Junction_info(k,13) > 0
            % get junction info
            ref_idx = sid2cid( Junction_info(k,2) );
            tgt_idx = sid2cid( Junction_info(k,3) );
            jlen = Junction_info(k,5);
            ref_pos = Junction_info(k,6);
            tgt_pos = Junction_info(k,9);
            sr = junction_matrix(ref_idx).mat(tgt_idx,ref_pos:1:ref_pos+jlen-1);
            st = junction_matrix(tgt_idx).mat(ref_idx,tgt_pos:1:tgt_pos+jlen-1);

            msk = sign( sr + st );
            % init junction_mat
            for m = 0:1:jlen-1
                if msk(m+1) == 0
                    junction_matrix(ref_idx).mat(tgt_idx,ref_pos+m) = tgt_pos+m;
                    junction_matrix(tgt_idx).mat(ref_idx,tgt_pos+m) = ref_pos+m;
                else
                    % junction_mat(tgt_idx,ref_pos+m,ref_idx) = 0;
                    % junction_mat(ref_idx,tgt_pos+m,tgt_idx) = 0;
                end
            end
        end
    end
    
    % fill junction_mat by message passing and check conflict
    n_loop_out = 0;
    while(1)
        % Try to fill companion
        n_loop = 0;
        while(1)
            n_conflict = 0;
            b_tmp = 0;
            for k = 1:1:n_seg_in_group
                for m = 1:1:seg(seg_idxs(k)).len
                    jvec = junction_matrix(k).mat(:,m);
                    for m1 = 1:1:n_seg_in_group
                        if m1 ~= k && jvec(m1) ~= 0
                            jvec2 = junction_matrix(m1).mat(:,jvec(m1));
                            for m2 = 1:1:n_seg_in_group
                                if m2 == k || m2 == m1

                                else
                                    if jvec(m2) == 0
                                        if jvec2(m2) ~= 0
                                            jvec(m2) = jvec2(m2);
                                            b_tmp = 1;
                                        end
                                    else
                                        if jvec2(m2) == 0
                                            jvec2(m2) = jvec(m2);
                                            b_tmp = 1;
                                        else
                                            if jvec(m2) == jvec2(m2)
                                                % OK
                                            else
                                                n_conflict = n_conflict + 1;
                                            end
                                        end
                                    end
                                end
                            end
                            junction_matrix(m1).mat(:,jvec(m1)) = jvec2;
                        end
                    end
                    junction_matrix(k).mat(:,m) = jvec;
                end
            end
            if b_tmp == 0
                break;
            else
                n_loop = n_loop +1;
            end
        end
        
        % If there were conflicts, resolve them
        if n_conflict > 0
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Try to Resolve conflict
            for k1 = 1:1:2
                n_conflict = 0;
                % b_tmp = 0;
                for k = 1:1:n_seg_in_group
                    for m = 1:1:seg(seg_idxs(k)).len
                        jvec = junction_matrix(k).mat(:,m);
                        for m1 = 1:1:n_seg_in_group
                            if m1 ~= k && jvec(m1) ~= 0
                                jvec2 = junction_matrix(m1).mat(:,jvec(m1));
                                for m2 = 1:1:n_seg_in_group
                                    if m2 == k || m2 == m1
                                    else
                                        if jvec(m2) == 0
                                            if jvec2(m2) ~= 0
                                                junction_matrix(m1).mat(m2,jvec(m1)) = 0;
                                            end
                                        else
                                            if jvec2(m2) == 0
                                                junction_matrix(k).mat(m2,m) = 0;
                                            else
                                                if jvec(m2) == jvec2(m2)
                                                    % OK
                                                else
                                                    % conflict
                                                    n_conflict = n_conflict + 1;
                                                    % resolve conflict
                                                    junction_matrix(m2).mat(k,jvec2(m2)) = 0;
                                                    junction_matrix(m2).mat(m1,jvec(m2)) = 0;
                                                    junction_matrix(k).mat(m1,m) = 0;
                                                    % junction_mat(m2,m,k) = 0;
                                                    junction_matrix(m1).mat(k,jvec(m1)) = 0;
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Check conflict again
            n_conflict = 0;
            % b_tmp = 0;
            for k = 1:1:n_seg_in_group
                for m = 1:1:seg(seg_idxs(k)).len
                    jvec = junction_matrix(k).mat(:,m);
                    for m1 = 1:1:n_seg_in_group
                        if m1 ~= k && jvec(m1) ~= 0
                            jvec2 = junction_matrix(m1).mat(:,jvec(m1));
                            for m2 = 1:1:n_seg_in_group
                                if m2 == k || m2 == m1
                                else
                                    if jvec(m2) == 0
                                        if jvec2(m2) ~= 0
                                            n_conflict = n_conflict + 1;
                                            junction_matrix(m1).mat(m2,jvec(m1)) = 0;
                                        end
                                    else
                                        if jvec2(m2) == 0
                                            n_conflict = n_conflict + 1;
                                            junction_matrix(k).mat(m2,m) = 0;
                                        else
                                            if jvec(m2) == jvec2(m2)
                                                % OK
                                            else
                                                % conflict
                                                n_conflict = n_conflict + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        n_loop_out = n_loop_out + 1;
        if n_conflict > 0
            str = sprintf('   WARNING: Conflicts found in Group %d ', group_idx );
            fprintf('\n   %s', str);
            fprintf(fp_log, '\n   %s', str);
            Nchar = 0;
        else
           break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set seg_vector using junction matrix %
    seg_vector = uint32( zeros(n_seg_in_group, max_seg_length) );
    n_seg_found_in_group = 0;
    for k = 1:1:n_seg_in_group
        cidx = k;
        sidx = seg_idxs(cidx);
        slen = seg(sidx).len;
        
        j_vec_old = zeros(n_seg_in_group,1); 
        for m = 1:1:slen
            
            j_vec_new = junction_matrix(cidx).mat(:,m);
            % Check if boundary
            b_bndry = 0;
            for m1 = 1:1:n_seg_in_group
                if j_vec_new(m1) == 0 
                    if j_vec_old(m1) == 0
                    else
                        b_bndry = 1;
                    end
                else
                    if j_vec_old(m1) == 0
                        b_bndry = 1;
                    else
                        if j_vec_new(m1) == (j_vec_old(m1)+1)

                        else
                            b_bndry = 1;
                        end
                    end
                end
            end
            if m == 1
                b_bndry = 1;
            end
            
            if seg_vector(k,m) == 0
                
                if b_bndry == 0
                    % check if there is any companion
                    % if not
                    if sum( abs(j_vec_new) ) == 0
                        seg_vector(k,m) = seg_vector(k,m-1);
                    else
                        % Get companion if any
                        n_tmp = 0;
                        for m1 = 1:1:length(j_vec_new)
                            if j_vec_new(m1) ~= 0
                                n_tmp = n_tmp + 1;
                            end
                        end
                        t_idxs = zeros(n_tmp,1);
                        n_tmp = 0;
                        for m1 = 1:1:length(j_vec_new)
                            if j_vec_new(m1) ~= 0
                                n_tmp = n_tmp + 1;
                                t_idxs(n_tmp) = m1;
                            end
                        end
                        % check if any companion is already written
                        wval = 0;
                        for m1 = 1:1:length(t_idxs)
                            if seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) == 0
                            else
                                if wval ~= 0 && wval ~= seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) )
                                    % fprintf('\n      ERROR: if wval ~= 0 && wval ~= seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) )');
                                end
                                wval = seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) );
                            end
                        end
                        if wval == 0
                            seg_vector(k,m) = seg_vector(k,m-1);
                            % set seg_vector for companion
                            for m1 = 1:1:length(t_idxs)
                                seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) = seg_vector(k,m);
                            end
                        else
                            seg_vector(k,m) = wval;
                            % set seg_vector for companion
                            for m1 = 1:1:length(t_idxs)
                                seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) = seg_vector(k,m);
                            end
                        end
                    end
                else % b_bndry > 0
                    % check if there is any companion
                    % if not
                    if sum( abs(j_vec_new) ) == 0
                        n_seg_found_in_group = n_seg_found_in_group + 1;
                        seg_vector(k,m) = n_seg_found_in_group;
                    else
                        % Get companion if any
                        n_tmp = 0;
                        for m1 = 1:1:length(j_vec_new)
                            if j_vec_new(m1) ~= 0
                                n_tmp = n_tmp + 1;
                            end
                        end
                        t_idxs = zeros(n_tmp,1);
                        n_tmp = 0;
                        for m1 = 1:1:length(j_vec_new)
                            if j_vec_new(m1) ~= 0
                                n_tmp = n_tmp + 1;
                                t_idxs(n_tmp) = m1;
                            end
                        end
                        % check if any companion is already written
                        wval = 0;
                        for m1 = 1:1:length(t_idxs)
                            if seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) == 0
                            else
                                if wval ~= 0 && wval ~= seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) )
                                    %  fprintf('\n      ERROR: if wval ~= 0 && wval ~= seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) )');
                                end
                                wval = seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) );
                            end
                        end
                        if wval == 0
                            n_seg_found_in_group = n_seg_found_in_group + 1;
                            seg_vector(k,m) = n_seg_found_in_group;
                            % set seg_vector for companion
                            for m1 = 1:1:length(t_idxs)
                                seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) = seg_vector(k,m);
                            end
                        else
                            seg_vector(k,m) = wval;
                            % set seg_vector for companion
                            for m1 = 1:1:length(t_idxs)
                                seg_vector( t_idxs(m1), j_vec_new(t_idxs(m1)) ) = seg_vector(k,m);
                            end
                        end
                    end
                end
                
            else % seg_vector(k,m) ~= 0
            end
            j_vec_old = j_vec_new;
        end % for m
    end % for k
    clear junction_matrix;
    
    b_segvec_valid = zeros(n_seg_found_in_group,1);
    for m2 = 1:1:n_seg_in_group
        s_tmp = seg_idxs(m2);
        l_tmp = seg(s_tmp).len;
        for m3 = 1:1:l_tmp
            b_segvec_valid( seg_vector(m2,m3) ) = 1;
        end
    end
    if sum(b_segvec_valid) ~= n_seg_found_in_group
        n_seg_found_in_group_tmp = 0;
        si_mapping = zeros(n_seg_found_in_group,1);
        for m2 = 1:1:n_seg_found_in_group
            if b_segvec_valid(m2) > 0
                n_seg_found_in_group_tmp = n_seg_found_in_group_tmp + 1;
                si_mapping(m2) = n_seg_found_in_group_tmp;
            end
        end
        for m2 = 1:1:n_seg_in_group
            s_tmp = seg_idxs(m2);
            l_tmp = seg(s_tmp).len;
            for m3 = 1:1:l_tmp
                seg_vector(m2,m3) = si_mapping( seg_vector(m2,m3) );
            end
        end
        n_seg_found_in_group = sum(b_segvec_valid);
    end
    
    seg_vector_all = seg_vector;
    for k = 1:1:n_seg_in_group
        sidx = seg_idxs(k);
        slen = seg(sidx).len;
        seg_vector(k, 1:slen) = seg_vector(k, 1:slen) + n_seg_found;
    end
    Group_Size_new(group_idx) = n_seg_found_in_group;
    
    % update n_seg_found_in_group
    n_seg_found = n_seg_found + n_seg_found_in_group;
    n_seg_added_in_group = n_seg_found_in_group - n_seg_in_group;
    n_seg_added = n_seg_added + n_seg_added_in_group;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create new segments and seg_connection_matrix
    if n_seg_found > Seg_mem_size
        n_seg_inc = 3000;
        seg_n = [seg_n; repmat( sg, n_seg_inc, 1 ) ];
        seg_connection_mat_lst = [seg_connection_mat_lst; uint32( zeros(n_seg_inc,3) )]; 
        group_id_new = [group_id_new; zeros(n_seg_inc,1)];
        Seg_mem_size = Seg_mem_size + n_seg_inc;
    end

    % get seg_vector for this group
    % seg_vector = seg_vector_all(seg_idxs, 1:max_seg_length);
    % Build seg_connection_mat and check loop
    seg_connection_mat = zeros( n_seg_found_in_group, n_seg_found_in_group, 'int8' );
    for k = 1:1:n_seg_in_group
        cidx = k;
        sidx = seg_idxs(cidx);
        slen = seg(sidx).len;

        idx_all = seg_vector_all(k,1);
        idx = seg_vector(k,1);
        group_id_new(idx) = group_idx;
        % pos = 1;
        len = 1;
        for m = 2:1:slen-1
            if seg_vector_all(k,m) == idx_all
                len = len + 1;
            else
                % Update seg lib and seg connection map
                seg_connection_mat(idx_all, seg_vector_all(k,m)) = 1;
                % Init variables
                idx_all = seg_vector_all(k,m);
                idx = seg_vector(k,m);
                group_id_new(idx) = group_idx;
                % pos = m;
                len = 1;
            end
        end % for m
        m = slen;
        if seg_vector_all(k,m) == idx_all
        else
            % Update seg lib and seg connection map
            seg_connection_mat(idx_all, seg_vector_all(k,m)) = 1;
        end
    end % for k
    
    group_cnt_pos(n,1) = seg_connection_mat_cnt + 1;
    for m1 = 1:1:n_seg_found_in_group
        for m2 = 1:1:n_seg_found_in_group
            if seg_connection_mat(m1,m2) > 0
                seg_connection_mat_cnt = seg_connection_mat_cnt + 1;
                seg_connection_mat_lst( seg_connection_mat_cnt, : ) = ...
                    [m1+n_seg_found-n_seg_found_in_group m2+n_seg_found-n_seg_found_in_group group_idx];
            end
        end
    end
    group_cnt_pos(n,2) = seg_connection_mat_cnt;
    group_cnt_pos(n,3) = group_idx;
    
    % get seg_vector for this group
    % seg_vector = seg_vector_all(seg_idxs, 1:max_seg_length);
    % Partitioning segments to create new set of segments
    for k = 1:1:n_seg_in_group
        cidx = k;
        sidx = seg_idxs(cidx);
        slen = seg(sidx).len;

        idx = seg_vector(k,1);
        group_id_new(idx) = group_idx;
        pos = 1;
        len = 1;
        for m = 2:1:slen-1
            if seg_vector(k,m) == idx
                len = len + 1;
            else
                % Update seg lib and seg connection map
                if seg_n(idx).len == 0
                    seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep = seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep ))) );
                else
                    if seg_n(idx).len < len
                        cd_tmp = uint32( zeros(4,len) );
                        cd_tmp(:,1:seg_n(idx).len) = seg_n(idx).cvg_dep(:,1:seg_n(idx).len);
                        seg_n(idx).len = uint32(len);
                        seg_n(idx).cvg_dep = cd_tmp + seg(sidx).cvg_dep(:,pos:pos+len-1);
                        [axx, mxi] = max( seg_n(idx).cvg_dep );
                        seg_n(idx).seq = int8( mxi-1 );
                        seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                    else
                        if seg_n(idx).len > len
                            % fprintf('PE21');
                        end
                        % seg_n(idx).len = uint32(len);
                        seg_n(idx).cvg_dep(:,1:len) = seg_n(idx).cvg_dep(:,1:len) + seg(sidx).cvg_dep(:,pos:pos+len-1);
                        [axx, mxi] = max( seg_n(idx).cvg_dep );
                        seg_n(idx).seq = int8( mxi-1 );
                        seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                    end
                end                
                % Init variables
                idx = seg_vector(k,m);
                group_id_new(idx) = group_idx;
                pos = m;
                len = 1;
            end
        end % for m
        m = slen;
        if seg_vector(k,m) == idx
            len = len + 1;
            if seg_n(idx).len == 0
                seg_n(idx).len = uint32(len);
                seg_n(idx).cvg_dep = seg(sidx).cvg_dep(:,pos:pos+len-1);
                [axx, mxi] = max( seg_n(idx).cvg_dep );
                seg_n(idx).seq = int8( mxi-1 );
                seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
            else
                if seg_n(idx).len < len
                    cd_tmp = uint32( zeros(4,len) );
                    cd_tmp(:,1:seg_n(idx).len) = seg_n(idx).cvg_dep(:,1:seg_n(idx).len);
                    seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep = cd_tmp + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                else
                    if seg_n(idx).len > len
                        % fprintf('PE22');
                    end
                    % seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep(:,1:len) = seg_n(idx).cvg_dep(:,1:len) + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                end
            end                
        else
            % Update seg lib and seg connection map
            if seg_n(idx).len == 0
                seg_n(idx).len = uint32(len);
                seg_n(idx).cvg_dep = seg(sidx).cvg_dep(:,pos:pos+len-1);
                [axx, mxi] = max( seg_n(idx).cvg_dep );
                seg_n(idx).seq = int8( mxi-1 );
                seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
            else
                if seg_n(idx).len < len
                    cd_tmp = uint32( zeros(4,len) );
                    cd_tmp(:,1:seg_n(idx).len) = seg_n(idx).cvg_dep;
                    seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep = cd_tmp + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                else
                    if seg_n(idx).len > len
                        % fprintf('PE23');
                    end
                    % seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep(:,1:len) = seg_n(idx).cvg_dep(:,1:len) + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                end
            end                
            
            % Init variables
            idx = seg_vector(k,m);
            group_id_new(idx) = group_idx;
            pos = m;
            len = 1;
            
            if seg_n(idx).len == 0
                seg_n(idx).len = uint32(len);
                seg_n(idx).cvg_dep = seg(sidx).cvg_dep(:,pos:pos+len-1);
                [axx, mxi] = max( seg_n(idx).cvg_dep );
                seg_n(idx).seq = int8( mxi-1 );
                seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
            else
                if seg_n(idx).len < len
                    cd_tmp = uint32( zeros(4,len) );
                    cd_tmp(:,1:seg_n(idx).len) = seg_n(idx).cvg_dep;
                    seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep = cd_tmp + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                else
                    if seg_n(idx).len > len
                        % fprintf('PE24');
                    end
                    % seg_n(idx).len = uint32(len);
                    seg_n(idx).cvg_dep(:,1:len) = seg_n(idx).cvg_dep(:,1:len) + seg(sidx).cvg_dep(:,pos:pos+len-1);
                    [axx, mxi] = max( seg_n(idx).cvg_dep );
                    seg_n(idx).seq = int8( mxi-1 );
                    seg_n(idx).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx).cvg_dep )) ));
                end
            end                
        end
    end % for k
    else
        seg_idxs = find( group_id(1:n_seg) == group_idx );
        if ~isempty(seg_idxs)
            n_seg_found = n_seg_found + 1;
            idx = n_seg_found;
            Group_Size_new(group_idx) = 1;
            group_id_new(idx) = group_idx;
            seg_n(idx) = seg(seg_idxs(1));
            n_seg_in_group = 1;
        else
            n_seg_in_group = 0;
        end
        group_cnt_pos(n,1) = seg_connection_mat_cnt + 1;
        group_cnt_pos(n,2) = seg_connection_mat_cnt;
        group_cnt_pos(n,3) = group_idx;
    end
    
    if n_seg_in_group > 1
        fprintf(repmat('\b', 1, Nchar));
        Nchar = fprintf( '%s %d/%d ', str_hdr, n, n_seg ); 
    else
        if mod(n, 100) == 0
            % fprintf('.');
            fprintf(repmat('\b', 1, Nchar));
            Nchar = fprintf( '%s %d/%d ', str_hdr, n, n_seg ); 
        end
    end
end

fprintf(repmat('\b', 1, Nchar));
str_disp = sprintf('%s %d -> %d added -> %d segs', str_hdr, n_seg, n_seg_added, n_seg_found );
fprintf( '%s', str_disp );
fprintf( fp_log, '%s', str_disp );
n_seg = n_seg_found;

clear seg;

%% Connect singly connected segments

str_disp = sprintf('   Compacting ... N_seg: %d ', n_seg );
fprintf( '\n%s', str_disp );
fprintf( fp_log, '\n%s', str_disp );
while(1)
    b_seg_valid = ones(1,n_seg);
    n_step = round(N_groups/3);
    Group_Size_new = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id_new(k);
        Group_Size_new(m) = Group_Size_new(m) + 1;
    end
    n_seg_count = 0;
    n_seg_found = 0;
    seg_connection_mat_cnt = 0;
    for n = 1:1:N_groups
        % group_idx = group_idx_sorted(n); %n; % Must not be sorted
        % [n group_cnt_pos(n,:) n_seg_count n_seg_found seg_connection_mat_cnt Group_Size_new(n)]
        sp = group_cnt_pos(n,1);
        ep = group_cnt_pos(n,2);
        % if sp > 0 && ep > 0 
        group_idx = group_cnt_pos(n,3);
        lst_t = double( seg_connection_mat_lst(sp:ep,:) );
        n_seg_selected = Group_Size_new( group_idx );
        seg_connection_mat_sel = sparse( lst_t(:,1)-n_seg_count, ...
                                 lst_t(:,2)-n_seg_count, ...
                                 lst_t(:,3), ...
                                 n_seg_selected, n_seg_selected );
        seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
        sidxs = (1:1:n_seg_selected) + n_seg_count;
        for k = 1:1:n_seg_selected
            kdx_r = k;
            idx_r = kdx_r + n_seg_count;
            if b_seg_valid(idx_r) > 0
                n_out_deg = sum( (seg_connection_mat_selected(kdx_r,:)).*b_seg_valid(sidxs) );
                if n_out_deg == 1
                    kdx_o = find( (seg_connection_mat_selected(kdx_r,:)).*b_seg_valid(sidxs) );
                    idx_o = kdx_o + n_seg_count;
                    if b_seg_valid(idx_o) > 0
                        n_in_deg = sum( (seg_connection_mat_selected(:,kdx_o)).*(b_seg_valid(sidxs)') );
                        if n_in_deg == 1
                            seg_len_new = double(seg_n(idx_r).len + seg_n(idx_o).len);
                            seg_n(idx_r).cvg_dep = [seg_n(idx_r).cvg_dep seg_n(idx_o).cvg_dep];
                            seg_n(idx_r).len = uint32(seg_len_new);
                            [axx, mxi] = max( seg_n(idx_r).cvg_dep );
                            seg_n(idx_r).seq = int8( mxi-1 );
                            seg_n(idx_r).ave_cvg_dep = uint32( round(mean( sum( seg_n(idx_r).cvg_dep )) ));

                            seg_connection_mat_selected(kdx_r,:) = seg_connection_mat_selected(kdx_o,:);
                            seg_connection_mat_selected(kdx_o,:) = int8( 0 );

                            b_seg_valid(idx_o) = 0;
                            if b_inter_disp > 0
                                str_disp = sprintf('Seg %d combined with Seg %d', idx_o, idx_r );
                                fprintf('%s', str_disp );
                            end
                        end
                    end
                end
            end
        end
        
        seg_cnt = 0;
        for k = 1:n_seg_selected
            if b_seg_valid(sidxs(k)) > 0
                seg_cnt = seg_cnt + 1;
                group_id_new(seg_cnt+n_seg_found) = group_id_new(k+n_seg_count);
                seg_n(seg_cnt+n_seg_found) = seg_n(k+n_seg_count);
                seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
            else
            end
        end
        
        group_cnt_pos(n,1) = seg_connection_mat_cnt + 1;
        for m1 = 1:1:seg_cnt
            for m2 = 1:1:seg_cnt
                if seg_connection_mat_selected(m1,m2) > 0
                    seg_connection_mat_cnt = seg_connection_mat_cnt + 1;
                    seg_connection_mat_lst( seg_connection_mat_cnt, : ) = ...
                        [m1+n_seg_found m2+n_seg_found group_idx];
                end
            end
        end
        group_cnt_pos(n,2) = seg_connection_mat_cnt;
        group_cnt_pos(n,3) = group_idx;
        
        n_seg_count = n_seg_count + n_seg_selected;
        n_seg_found = n_seg_found + seg_cnt;
        % end
        if mod(n, n_step) == 0
            fprintf('>' );
            fprintf(fp_log, '>' );
        end
    end
    n_seg = n_seg_found;
    Group_Size_new = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = group_id_new(k);
        Group_Size_new(m) = Group_Size_new(m) + 1;
    end
    
    if sum( b_seg_valid(1:n_seg) ) ~= n_seg 
        str_disp = sprintf(' %d ', n_seg );
        fprintf( '%s', str_disp );
        fprintf( fp_log, '%s', str_disp );
    else
        str_disp = sprintf(' %d ', n_seg );
        fprintf( '%s', str_disp );
        fprintf( fp_log, '%s', str_disp );
        break;
    end
end

%% Save data
str_disp = sprintf('   Saving segment info. ' );
fprintf( '\n%s', str_disp );
fprintf( fp_log, '\n%s', str_disp );

Group_index = group_id_new(1:n_seg);

fname_txt = sprintf('%s.sgmnt', fname_out );
fp_t = fopen( fname_txt, 'wt' );

kmod = round(n_seg/8);
for k = 1:1:n_seg
    fprintf(fp_t,'>Segment\t%d\t%d\t%d\n', k, Group_index(k), seg_n(k).len );
    fprintf(fp_t,'%s\n', sub_NumSeq2NTstr( seg_n(k).seq(1:seg_n(k).len) ) );
    for m2 = 1:1:4
        for m1 = 1:1:seg_n(k).len-1
            fprintf(fp_t,'%d\t', seg_n(k).cvg_dep(m2,m1) );
        end
        m1 = seg_n(k).len;
        fprintf(fp_t,'%d\n', seg_n(k).cvg_dep(m2,m1) );
    end
            
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
fclose(fp_t);

fname_txt2 = sprintf('%s.scmat', fname_out );
fp_t = fopen( fname_txt2, 'wt' );
fprintf(fp_t,'%d\t%d\t%d\t%d\n', n_seg, N_groups, N_reads_valid, N_reads_total );
kmod = round(seg_connection_mat_cnt/8);
for k = 1:1:seg_connection_mat_cnt
    fprintf(fp_t,'%d\t%d\n', seg_connection_mat_lst(k,1), seg_connection_mat_lst(k,2) );
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
fclose(fp_t);

% fname_out2 = sprintf('%s', out_file_prefix );
% fname_txt = sprintf('%s.rcmap3', fname_out2 );
% fp_t = fopen( fname_txt, 'wt' );
% fprintf(fp_t, '%d\n', cntg_cnt);
% kmod = round(cntg_cnt/6);
% for k = 1:1:cntg_cnt
%     fprintf(fp_t, '%d\n', cntg_valid_ind(k));
%     if mod(k,kmod) == 0
%         fprintf('.');
%     end
% end
% Nrt = round( N_reads_total/(R_mode+1) );
% fprintf(fp_t, '%d\n', Nrt);
% kmod = round(Nrt/6);
% for k = 1:1:Nrt
%     fprintf(fp_t, '%d\t%d\n', read_cntg_map(k,1), read_cntg_map(k,2) );
%     if mod(k,kmod) == 0
%         fprintf('.');
%     end
% end
% fclose(fp_t);

str_disp = sprintf(' saved to %s/scmat', fname_txt );
fprintf( '%s\n', str_disp );
fprintf( fp_log, '%s\n', str_disp );

fprintf('%s', dstr);
str_disp = sprintf(' and completed at %s', datestr(now) );
fprintf('%s\n', str_disp);

fprintf(fp_log, '%s', dstr );
fprintf(fp_log, '  and completed %s\n\n', datestr(now) );
fclose(fp_log);

end

%% %%%%%%%%% %%
%% functions %%
%% %%%%%%%%% %%

function seq_int8 = f03_seq_est( seg_cvg )
    [axx, mxi] = max( seg_cvg );
    seq_int8 = int8( mxi-1 );
end

% function cvg = f03_set_cvg( seg_seq )
%     seg_len = length( seg_seq );
%     cvg = uint16( zeros(4,seg_len) );
%     for k = 1:1:seg_len
%         cvg(seg_seq(k)+1,k) = 1;
%     end
% end
% 
% function tvec = f04_check_loop( n_seg_selected, seg_connection_mat_selected )
%     % b_tmp = 0;
%     Tmat = eye(n_seg_selected);
%     for k = 1:1:40
%         Tmat = sign(Tmat*double( seg_connection_mat_selected(1:n_seg_selected,1:n_seg_selected)));
%         tvec = diag(Tmat); 
%         if sum( tvec ) > 0
%             % b_tmp = 1;
%             break;
%         end
%     end
% end
% 
%% f04_check_polyA_tail_Num_v03

function [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04a( seq, Min_tail_length, d_threshold, n_rd_len )

    persistent d_th;
    if isempty( d_th )
        d_th = ceil( (1:1:n_rd_len).*d_threshold );
    end
    % tail_min_len = Min_tail_length;
    Len = length(seq);
    %d_th = d_threshold;
    
    if Len > Min_tail_length
        ts = sum( abs(sign( seq(1:Min_tail_length) - 'T' )) );
        if ts > d_th(Min_tail_length);
            tail_len_b = 0;
        else
            b_tmp = 0;
            for k = 1:1:Len-Min_tail_length
                if seq(k+Min_tail_length) ~= 'T'
                    ts = ts + 1;
                    if ts > d_th(k+Min_tail_length) %ceil( (k+1)*pe )
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

        % Tseq = sign(seq);
        b_tmp = 0;
        ts = sum( abs(sign( seq(Len-Min_tail_length+1:Len) - 'A' )) );
        if ts > d_th(Min_tail_length)
            tail_len_f = 0;
        else
            for k = 1:1:Len-Min_tail_length
                if seq(Len-Min_tail_length+1-k) ~= 'A'
                    ts = ts + 1;
                    if ts > d_th(k+Min_tail_length) %ceil( (k+1)*pe )
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
            % tail_len = 0;
            tail_len_f = 0;
            tail_len_b = 0;
        else
            if ts1 <= d_th(Len)
                % tail_len = -Len;
                tail_len_f = Len;
                tail_len_b = Len;
            else
                % tail_len = Len;
                tail_len_f = Len;
                tail_len_b = Len;
            end
        end
    end
end


%% f04_strand_correction( sc_ref, mode )
function [dir, sc_k] = f04_strand_correction( sc_ref, mode )
    switch( sc_ref )
        case 1,
            switch( mode )
                case +1,  dir = -1; sc_k = +1;
                case +2,  dir = +1; sc_k = +1;
                case -1,  dir = -1; sc_k = -1;
                case -2,  dir = +1; sc_k = -1;
                otherwise,  dir = 0; sc_k = 0;
            end
        case -1,
            switch( mode )
                case +1,  dir = +1; sc_k = -1;
                case +2,  dir = -1; sc_k = -1;
                case -1,  dir = +1; sc_k = +1;
                case -2,  dir = -1; sc_k = +1;
                otherwise,  dir = 0; sc_k = 0;
            end
        otherwise, dir = 0; sc_k = 0;
    end
end

%% Compatibility check

function [b_compatibility, c_mode_lst] = f04_compatibility_check( junction_info_tmp1, junction_info_tmp2, b_inter_disp )

    junc_mode = zeros(2,1);
    junc_mode(1) = junction_info_tmp1(12);
    junc_mode(2) = junction_info_tmp2(12);
    % check compatibility
    b_compatibility = [0 0];
    c_mode_lst = [];
    switch( junc_mode(1) )
        case 1,
            switch( junc_mode(2) )
                case 1,  
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 21');
                    end
                case 2,     
                    if junction_info_tmp1(6) < junction_info_tmp2(6)
                        b_compatibility = [1 1];    
                        c_mode_lst = [2 1];
                    else
                        if junction_info_tmp1(5)-junction_info_tmp1(14) > junction_info_tmp2(5)-junction_info_tmp2(14)
                            b_compatibility = [1 0];    
                            c_mode_lst = [2 1];
                        else
                            b_compatibility = [0 1];    
                            c_mode_lst = [1 2];
                        end
                    end
                case 3,    
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 22');
                    end
                case 4,     
                    if b_inter_disp > 0
                    end
                    if junction_info_tmp1(5)-junction_info_tmp1(14) > junction_info_tmp2(5)-junction_info_tmp2(14)
                        b_compatibility = [1 0];
                        c_mode_lst = [2 1];
                    else
                        b_compatibility = [0 1];
                        c_mode_lst = [1 2];
                    end
                otherwise,  
                    b_compatibility = [1 0];
                    c_mode_lst = [2 1];
            end
        case 2, 
            if b_inter_disp > 0
                fprintf('\n      ERROR: This case cannot happen 11');
            end
        case 3,
            switch( junc_mode(2) )
                case 1,   
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 23');
                    end
                case 2,  
                    if b_inter_disp > 0
                    end
                    if junction_info_tmp1(5)-junction_info_tmp1(14) > junction_info_tmp2(5)-junction_info_tmp2(14)
                        b_compatibility = [1 0];
                        c_mode_lst = [-1 -1];
                    else
                        b_compatibility = [0 1];
                        c_mode_lst = [-2 -2];
                    end
                case 3,     
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 24');
                    end
                case 4,     
                    if junction_info_tmp1(6) > junction_info_tmp2(6)
                        b_compatibility = [1 1];    
                        c_mode_lst = [-2 -2];
                    else
                        if junction_info_tmp1(5)-junction_info_tmp1(14) > junction_info_tmp2(5)-junction_info_tmp2(14)
                            b_compatibility = [1 0];    
                            c_mode_lst = [-1 -1];
                        else
                            b_compatibility = [0 1];    
                            c_mode_lst = [-2 -2];
                        end
                    end
                otherwise,  
                    b_compatibility = [1 0];
                    c_mode_lst = [-1 -1];
            end
        case 4,
            if b_inter_disp > 0
                fprintf('\n      ERROR: This case cannot happen 12');
            end
        otherwise,
            switch( junc_mode(2) )
                case 1,   
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 25');
                    end
                case 2,     
                    b_compatibility = [0 1];    
                    c_mode_lst = [1 2];
                case 3,    
                    if b_inter_disp > 0
                        fprintf('\n      ERROR: This case cannot happen 26');
                    end
                case 4,     
                    b_compatibility = [0 1];    
                    c_mode_lst = [-2 -2];
               otherwise,  
                    b_compatibility = [0 0];
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
