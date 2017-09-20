
function fname_fasta = trareco_td( config_file )

dstr = sprintf('Isoform detection started %s', char( datestr(now) ) );
disp(dstr);
if exist('b_all', 'var') == 0
    b_all = 1;
end

%% Load data
if ischar(config_file)
    cfg = trareco_get_config( config_file );
    if isnumeric(cfg)
        return;
    end
else
    return;
end
out_file_prefix = cfg.output_prefix;

R_mode = cfg.read_mode;
if R_mode == 0
    input_1 = cfg.input_file_1; 
else
    input_1 = cfg.input_file_1; 
    input_2 = cfg.input_file_2;
end

%% Set parameters

Cvg_dep_wgt_div = round( cfg.nominal_read_length*2 );
Edge_seg_len_min = 0; 
E_cwgt = 0;

Cset_size_max = 1000;
N_Csets_max = cfg.max_num_csets;
n_max_paths = cfg.max_num_paths;
b_tail_suppress = cfg.b_tail_suppress;
mse_mf = (100+cfg.mse_margin_percent)/100;
Min_seg_len_to_combine = 120;
Junction_overlap_backoff = cfg.junction_backoff;
% Max_num_paths_to_search = cfg.max_num_isoforms;
b_disp = 0;

%% Load SC matrix
if ~isempty( cfg.output_dir ) && ischar( cfg.output_dir )
    if ~exist( cfg.output_dir, 'dir' )
        fprintf( '   WARNING: The directory you specified not exists. \n' );
        return;
    end
    out_file_prefix = sprintf('%s/%s', cfg.output_dir, out_file_prefix ); 
end
fname_txt2 = sprintf('%s.scmat', out_file_prefix );
fp_t = fopen( fname_txt2, 'rt' );
linea = fgets(fp_t);
data = sscanf(linea, '%d', 4 );
seg_connection_mat_lst_tmp = fscanf(fp_t, '%d', [2, inf] )';
fclose(fp_t);

N_seg = data(1);
N_groups = data(2);
N_reads_valid = data(3);
N_reads_total = data(4);
[seg_connection_mat_cnt, nc] = size(seg_connection_mat_lst_tmp);
seg_connection_mat_lst = [seg_connection_mat_lst_tmp ones(seg_connection_mat_cnt,1)];
% Output: N_seg, N_groups
% Output: seg_connection_mat_cnt, seg_connection_mat_lst

%% Load Segments
fname_stlb = sprintf('%s.scmat/sgmnt', out_file_prefix );
str_disp = sprintf('   Loading .. %s ', fname_stlb );
fprintf('%s', str_disp);

n_seg = N_seg; 
sg.len = uint32(0);
sg.seq = int8( [] );
sg.cvg_dep = uint32( zeros(4,0) );
sg.group_index = uint32(0);
sg.ave_cvg_dep = 0;
seg = repmat( sg, n_seg, 1 );

fname_txt = sprintf('%s.sgmnt', out_file_prefix );
fp_t = fopen( fname_txt, 'rt' );
kmod = round(n_seg/10);
n_cnt = 0;
for k = 1:1:n_seg    
    linea = fgets(fp_t);
    if linea < 0
        break;
    else
        if linea(1) == '>'
            n_cnt = n_cnt + 1;
            
            ptr = 1;
            [a, b, c, next_idx] = sscanf( linea(ptr:end), '%s', 1 );
            ptr = ptr + next_idx;
            [val, a, b, c] = sscanf( linea(ptr:end), '%d', 3 );
            seg(n_cnt).len = val(3);
            seg(n_cnt).group_index = val(2);
            
            lineb = fgets(fp_t);
            if seg(n_cnt).len ~= (length(lineb)-1)
                fprintf('\n      ERROR occurred when reading file (%d) %d ~= %d ', n_cnt, seg(n_cnt).len, (length(lineb)-1));
            else
                seg(n_cnt).len = length(lineb)-1;
            end
            seg(n_cnt).seq = sub_NTstr2NumSeq( lineb(1:end-1) );
            seg(n_cnt).cvg_dep = ( zeros( 4, seg(n_cnt).len, 'uint32' ) );
            for m = 1:1:4
                linec = fgets(fp_t);
                consen_tmp = sscanf( linec, '%d', [1 seg(n_cnt).len] );
                seg(n_cnt).cvg_dep(m,:) = consen_tmp;
            end
            seg(n_cnt).ave_cvg_dep = mean( sum( seg(n_cnt).cvg_dep ) );
        end
    end
    
    if mod(k,kmod) == 0
        fprintf('.');
    end
end
fclose(fp_t);
str_disp = sprintf(' %d segments read', n_seg );
fprintf( '%s', str_disp );
% Output: seg(len, seq, cvg_dep, ave_cvg_dep, group_index), n_seg

Group_Size = zeros(N_groups,1);
for k = 1:1:n_seg
    m = seg(k).group_index;
    Group_Size(m) = Group_Size(m) + 1;
end
max_gsize = max( Group_Size );
Group_Size = zeros(N_groups,1);
Group_members = zeros( N_groups, max_gsize );
for k = 1:1:n_seg
    m = seg(k).group_index;
    Group_Size(m) = Group_Size(m) + 1;
    Group_members(m,Group_Size(m)) = k;
end

%% Log file
if b_all == 0
    % fname_tr = sprintf('%s_td%s_all', out_file_prefix, tver );
    fname_tr = sprintf('%s_all', out_file_prefix );
else
    % fname_tr = sprintf('%s_td%s_sel', out_file_prefix, tver ); %, Cvg_dep_wgt_div );
    fname_tr = sprintf('%s', out_file_prefix ); %, Cvg_dep_wgt_div );
end

fname_log = sprintf('%s.log', out_file_prefix );
fp_log = fopen( fname_log, 'a' );

if fp_log >= 0
    fprintf(fp_log, '%s\n', dstr);
    fprintf(fp_log, '   # Input: %s.sgmnt(scmat)', out_file_prefix );
else
    if cfg.b_split_merged == 0
    %     fprintf(fp_log, '(Possibly) Merged Path split' );
    else
    end
    fprintf(fp_log, '\n%s\n', dstr);
    fprintf(fp_log, '   # Input: %s.sgmnt(scmat)', out_file_prefix );
end
seg_cnt = 0;
N_bases = 0;
for n = 1:n_seg
    seg_cnt = seg_cnt + 1;
    seg(seg_cnt).len = seg(n).len;
    seg(seg_cnt).cvg_dep = seg(n).cvg_dep;
    seg(seg_cnt).group_index = seg(n).group_index;
    seg(seg_cnt).seq = f03_seq_est( seg(seg_cnt).cvg_dep );
    n_bases = sum( sum( seg(seg_cnt).cvg_dep ) );
    N_bases = N_bases + n_bases;
    seg(seg_cnt).ave_cvg_dep = n_bases/seg(seg_cnt).len;
end
n_seg = seg_cnt;
str_disp = sprintf(' .. %d bases', N_bases );
fprintf( '%s', str_disp );

%% Check Hidden path
% if cfg.b_split_merged > 0
%     str_disp = sprintf('   Checking merged path ' );
%     fprintf( '\n%s', str_disp );
% 
%     n_ave = zeros(1,n_seg);
%     n_stdev = zeros(1,n_seg);
%     p_ave = zeros(1,n_seg);
%     error_rate = zeros(1,n_seg);
%     e_threshold = zeros(1,n_seg);
%     n_split = 0;
%     n_seg_add = 0;
%     for k = 1:1:n_seg
%         idx_r = k;
%         slen = seg(idx_r).len;
%         cd_tmp = double( seg(idx_r).cvg_dep(:,1:slen) );
%         [mxv, mxi] = max( cd_tmp );
%         for m = 1:1:slen
%             cd_tmp(mxi(m),m) = 0;
%         end
%         n_ave(k) = mean( mean(cd_tmp) )*(4/3);
%         n_stdev(k) = sqrt( var( reshape(cd_tmp,1,[]) ) );
%         p_ave(k) = mean( mxv );
%         error_rate(k) = n_ave(k)/p_ave(k);
%         e_threshold(k) = max( min( p_ave(k)/8, ceil(n_ave(k))*6 ), cfg.min_cvg_depth_js );
% 
%         s_vector = zeros(1,slen);
%         for m = 1:1:slen
%             tv = max( cd_tmp(:,m) - e_threshold(k), 0 );
%             if sum(tv) > 0
%                 % [k m cvg_dep_lib_n(:,m,idx_r)' tv']
%                 [mv, s_vector(m)] = max(tv);
%                 n_split = n_split + 1;
%             end
%         end
%         if sum(s_vector) > 0
%             a_vector = sign(s_vector);
%             a_prev = a_vector(slen);
%             m_prev = slen;
%             for m = slen-1:-1:1
%                 if a_vector(m) ~= a_prev
%                     if a_prev == 0
%                         n_seg_add = n_seg_add + 1;
%                         idx = n_seg + n_seg_add;
%                         
%                         seg(idx).len = int16(m_prev - m);
%                         seg(idx).cvg_dep = seg(idx_r).cvg_dep(:,(m+1):m_prev);
%                         seg(idx).seq = f03_seq_est( seg(idx).cvg_dep );
%                         seg(idx).ave_cvg_dep = int16( round(mean( sum( seg(idx).cvg_dep )) ));
%                         
%                         group_id_new(idx) = group_id_new(idx_r);
% 
%                         seg(idx_r).len = int16(m);
%                         seg(idx_r).cvg_dep = seg(idx_r).cvg_dep(:,1:m);
%                         seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
%                         seg(idx_r).ave_cvg_dep = int16( round(mean( sum( seg(idx_r).cvg_dep )) ));
% 
%                         seg_connection_mat(1:idx,idx) = 0;
%                         seg_connection_mat(idx,1:idx) = 0;
%                         seg_connection_mat(idx,1:idx) = seg_connection_mat(idx_r,1:idx);
%                         seg_connection_mat(idx_r,1:idx) = 0;
%                         seg_connection_mat(idx_r,idx) = 1;
%                     else
%                         n_seg_add = n_seg_add + 1;
%                         idx = n_seg + n_seg_add;
%                         
%                         seg(idx).len = int16(m_prev - m);
%                         [mxv, mxi] = max( seg(idx_r).cvg_dep(:,(m+1):m_prev) );
%                         seg(idx).cvg_dep = ( zeros(4,seg(idx).len, 'uint32') );
%                         for m1 = 1:1:seg(idx).len
%                             seg(idx).cvg_dep(mxi,m1) = mxv; 
%                         end
%                         seg(idx).seq = f03_seq_est( seg(idx).cvg_dep );
%                         seg(idx).ave_cvg_dep = int16( round(mean( sum( seg(idx).cvg_dep )) ));
%                         group_id_new(idx) = group_id_new(idx_r);
%                         idx1 = idx;
% 
%                         n_seg_add = n_seg_add + 1;
%                         idx = n_seg + n_seg_add;
%                         seg(idx).len = int16(m_prev - m);
%                         [mxv, mxi] = max( seg(idx_r).cvg_dep(:,(m+1):m_prev) );
%                         seg(idx).cvg_dep = seg(idx_r).cvg_dep(:,(m+1):m_prev);
%                         for m1 = 1:1:seg_len_n(idx)
%                             seg(idx).cvg_dep(mxi,m1) = 0; 
%                         end
%                         seg(idx).seq = f03_seq_est( seg(idx).cvg_dep );
%                         seg(idx).ave_cvg_dep = int16( round(mean( sum( seg(idx).cvg_dep )) ));
%                         group_id_new(idx) = group_id_new(idx_r);
%                         idx2 = idx;
% 
%                         seg(idx_r).len = int16(m);
%                         seg(idx_r).cvg_dep = seg(idx_r).cvg_dep(:,1:seg(idx_r).len);
%                         seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
%                         seg(idx_r).ave_cvg_dep = int16( round(mean( sum( seg(idx_r).cvg_dep )) ));
%                         
%                         seg_connection_mat(1:idx2,idx1) = 0;
%                         seg_connection_mat(idx1,1:idx2) = 0;
%                         seg_connection_mat(1:idx2,idx2) = 0;
%                         seg_connection_mat(idx2,1:idx2) = 0;
%                         seg_connection_mat(idx1,1:idx2) = seg_connection_mat(idx_r,1:idx2);
%                         seg_connection_mat(idx2,1:idx2) = seg_connection_mat(idx_r,1:idx2);
%                         seg_connection_mat(idx_r,1:idx2) = 0;
%                         seg_connection_mat(idx_r,idx1) = 1;
%                         seg_connection_mat(idx_r,idx2) = 1;
%                     end
%                     a_prev = a_vector(m);
%                     m_prev = m;
%                 end
%             end
%         end
%         n_step = round(n_seg/5);
%         if mod(k, n_step) == 0
%             fprintf('.' );
%         end
%     end
%     str_disp = sprintf(' N_splt: %d, N_seg: %d -> %d ', n_split, n_seg, n_seg + n_seg_add );
%     fprintf( '%s\n', str_disp );
%     str_disp = sprintf('Na: %f, Pa: %f -> Estimated Error Rate: %f ', mean(n_ave), mean(p_ave), mean(error_rate) );
%     fprintf( '%s', str_disp );
%     n_seg = n_seg + n_seg_add;
% end

%% Suppressing multiple tail sites
if b_tail_suppress > 0
    str_disp = sprintf('   Removing tail ..... ' );
    fprintf( '\n%s', str_disp );
    fprintf( fp_log, '\n%s', str_disp );
    str_disp = sprintf('%d', n_seg );
    fprintf( '%s', str_disp );
    fprintf( fp_log, '%s', str_disp );
    Nchar = 0;
                
    b_seg_valid = ones(1,n_seg);
    % n_step = round(N_groups/3);
    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
    end
    max_gsize = max( Group_Size );
    Group_Size = zeros(N_groups,1);
    Group_members = zeros( N_groups, max_gsize );
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
        Group_members(m,Group_Size(m)) = k;
    end
    
    scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
    seg_n = repmat( sg, n_seg, 1 );
    seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
    seg_connection_mat_cnt_n = 0;
    n_seg_found = 0;
    n_seg_count = 0;
    for n = 1:1:N_groups
        % group_idx = group_idx_sorted(n); % Must not be sorted
        group_idx = n;
        if Group_Size(group_idx) == 1
            sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
            if length(sidxs) ~= 1
                fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
            end
            seg_n(1+n_seg_found) = seg(sidxs(1));
            n_seg_found = n_seg_found + 1;
            n_seg_count = n_seg_count + 1;
        else
            lst_sel = find( scm_gidx == group_idx );
            sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
            n_seg_selected = length(sidxs);
            if isempty(lst_sel)
                seg_connection_mat_selected = zeros(n_seg_selected);
            else
                tmp_idxs = zeros(n_seg,1);
                tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                seg_connection_mat_selected = int8( zeros(n_seg_selected) );
                for k = 1:1:length(lst_sel)
                    seg_connection_mat_selected( tmp_idxs(seg_connection_mat_lst(lst_sel(k),1)), ...
                        tmp_idxs(seg_connection_mat_lst(lst_sel(k),2)) ) = 1;
                end
            end
            if n_seg_selected ~= length(sidxs)
                fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
            end

            for k = 1:1:n_seg_selected
                idx_r = sidxs(k); 
                if b_seg_valid(idx_r) > 0 && seg(idx_r).len > Junction_overlap_backoff && seg(idx_r).len <= Min_seg_len_to_combine
                    n_in_deg = sum( double( ( seg_connection_mat_selected(:,k))).*b_seg_valid(sidxs)' );
                    n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                    if n_out_deg == 0 
                        if n_in_deg > 0
                            if sum( seg(idx_r).seq(Junction_overlap_backoff+1:seg(idx_r).len) ) == 0
                                idxs = find( double( ( seg_connection_mat_selected(:,k))) );
                                for m = 1:1:length(idxs)
                                    seg_connection_mat_selected(idxs(m),k) = 0;
                                end
                                b_seg_valid(idx_r) = 0;
                            end
                        end
                    else
                        if n_in_deg == 0
                            if sum( abs( seg(idx_r).seq(1:seg(idx_r).len-Junction_overlap_backoff) - 3 ) ) == 0
                                idxs = find( double( ( seg_connection_mat_selected(k,:))) );
                                for m = 1:1:length(idxs)
                                    seg_connection_mat_selected(k,idxs(m)) = 0;
                                end
                                b_seg_valid(idx_r) = 0;
                            end
                        end
                    end
                end
            end

            seg_cnt = 0;
            for k = 1:n_seg_selected
                if b_seg_valid(sidxs(k)) > 0
                    seg_cnt = seg_cnt + 1;
                    seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                    seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                    seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                else
                end
            end
            for m1 = 1:1:seg_cnt
                for m2 = 1:1:seg_cnt
                    if seg_connection_mat_selected(m1,m2) > 0
                        seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                        seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                            [m1+n_seg_found m2+n_seg_found group_idx];
                    end
                end
            end
            n_seg_found = n_seg_found + seg_cnt;
            n_seg_count = n_seg_count + n_seg_selected;
        end
        if mod( n, 20 ) == 0 || n == N_groups 
        if Nchar > 0
            fprintf(repmat('\b', 1, Nchar));
        end
        Nchar = fprintf( '/%d', n_seg_count );
        end
    end
    n_seg = n_seg_found;
    seg = seg_n(1:n_seg);
    seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
    seg_connection_mat_cnt = seg_connection_mat_cnt_n;

    Group_Size = zeros(N_groups,1);
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
    end
    max_gsize = max( Group_Size );
    Group_Size = zeros(N_groups,1);
    Group_members = zeros( N_groups, max_gsize );
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
        Group_members(m,Group_Size(m)) = k;
    end
    
    if sum( b_seg_valid(1:n_seg) ) ~= n_seg 
        if Nchar > 0
            fprintf(repmat('\b', 1, Nchar));
            Nchar = 0;
        end
        str_disp = sprintf(' -> %d', n_seg );
        fprintf( '%s', str_disp );
        fprintf( fp_log, '%s', str_disp );
    else
    end

    %% Connect singly connected segments
    % str_disp = sprintf('   Condensing ... %d ', n_seg );
    % fprintf( '\n%s', str_disp );
    n_loop = 1;
    while(1)
        b_seg_valid = ones(1,n_seg);
        % n_step = round(N_groups/3);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) == 1
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                if length(sidxs) ~= 1
                    fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                end
                seg_n(1+n_seg_found) = seg(sidxs(1));
                n_seg_found = n_seg_found + 1;
                n_seg_count = n_seg_count + 1;
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end

                for k = 1:1:n_seg_selected
                    idx_r = sidxs(k); 
                    if b_seg_valid(idx_r) > 0
                        n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                        if n_out_deg == 1
                            idx_o = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            if b_seg_valid(sidxs(idx_o)) > 0
                                n_in_deg = sum( double( ( seg_connection_mat_selected(:,idx_o))).*b_seg_valid(sidxs)' );
                                if n_in_deg == 1
                                    seg_len_new = double(seg(idx_r).len + seg(sidxs(idx_o)).len);
                                    seg(idx_r).cvg_dep = [seg(idx_r).cvg_dep seg(sidxs(idx_o)).cvg_dep];
                                    seg(idx_r).len = seg_len_new;
                                    [v, mxi] = max( seg(idx_r).cvg_dep );
                                    seg(idx_r).seq = int8( mxi-1 );
                                    seg(idx_r).ave_cvg_dep = ( ( mean( sum( seg(idx_r).cvg_dep ) )));

                                    seg_connection_mat_selected(k,:) = seg_connection_mat_selected(idx_o,:);
                                    seg_connection_mat_selected(idx_o,:) = 0; 

                                    b_seg_valid(sidxs(idx_o)) = 0;
                                end
                            end
                        end
                    end
                end
                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found group_idx];
                        end
                    end
                end
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            if mod( n, 20 ) == 0 || n == N_groups
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf( '/%d', n_seg_count );
            end
        end
        n_seg = n_seg_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;

        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
        
        if sum( b_seg_valid ) ~= length( b_seg_valid ) 
            n_loop = n_loop + 1;
%             str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%             fprintf( '%s', str_disp );
%             fprintf( fp_log, '%s', str_disp );
        else
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            str_disp = sprintf(' -> %d', n_seg );
            fprintf( '%s', str_disp );
            fprintf( fp_log, '%s', str_disp );
            break;
        end
    end
end

if b_tail_suppress > 0
    % str_disp = sprintf('   Suppressing multiple tail sites ...  N_seg: %d', n_seg );
    % fprintf( '\n%s', str_disp );
    for kk = 1:1:3
        b_seg_valid = ones(1,n_seg);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) == 1
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                if length(sidxs) ~= 1
                    fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                end
                seg_n(1+n_seg_found) = seg(sidxs(1));
                n_seg_found = n_seg_found + 1;
                n_seg_count = n_seg_count + 1;
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end

                for k = 1:1:n_seg_selected
                    idx_r = sidxs(k); 
                    if (b_seg_valid(idx_r) > 0) && (seg(idx_r).len <= Min_seg_len_to_combine)
                        n_in_deg = sum( double( ( seg_connection_mat_selected(:,k))).*b_seg_valid(sidxs)' );
                        n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                        if (n_out_deg == 0 && n_in_deg == 1) %|| (n_in_deg == 0 && n_out_deg == 1) 
                            idx_parent = find( double( ( seg_connection_mat_selected(:,k))).*b_seg_valid(sidxs)' );
                            idx_children = find( double( ( seg_connection_mat_selected(idx_parent(1),:))).*b_seg_valid(sidxs) );
                            n_children = length(idx_children);
                            % check if it is one stage tree
                            b_tmp = 0;
                            for m = 1:1:n_children
                                if seg( sidxs(idx_children(m)) ).len > Min_seg_len_to_combine || b_seg_valid(sidxs(idx_children(m))) == 0
                                    b_tmp = 1;
                                end
                                if sum( double( ( seg_connection_mat_selected(:,idx_children(m)))).*b_seg_valid(sidxs)' ) > 1
                                    b_tmp = 1;
                                else
                                    if sum( double( ( seg_connection_mat_selected(idx_children(m),:))).*b_seg_valid(sidxs) ) > 0
                                        b_tmp = 1;
                                    end
                                end
                            end
                            if b_tmp == 0 && n_children > 1
                                b_tmp = 0;
                                slen_eff = zeros(n_children,1);
                                for m = 1:1:n_children
                                    idx = sidxs(idx_children(m)); 
                                    [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04b( sub_NumSeq2NTstr( seg(idx).seq ), 0, 0 );
                                    if tail_len_b > 0
                                        b_tmp = 1;
                                    else
                                        slen_eff(m) = seg(idx).len - tail_len_f;
                                    end
                                end
                                if b_tmp == 0
                                    [v, mx_idx] = max( slen_eff );
                                    idx_l = sidxs(idx_children(mx_idx)); 
                                    for m = 1:1:n_children
                                        if m ~= mx_idx
                                            idx_s = sidxs(idx_children(m)); 
                                            if sum( abs(seg(idx_l).seq(1:slen_eff(m)) - seg(idx_s).seq(1:slen_eff(m))) ) > 0
                                                b_tmp = 1;
                                            end
                                        end
                                    end
                                    if b_tmp == 0
                                        for m = 1:1:n_children
                                            if m ~= mx_idx
                                                idx_s = sidxs(idx_children(m)); 
                                                sp = 1;
                                                ep = slen_eff(m); 
                                                seg(idx_l).cvg_dep(:,sp:ep) = ...
                                                   seg(idx_l).cvg_dep(:,sp:ep) + seg(idx_s).cvg_dep(:,sp:ep);
                                                [v, mxi] = max( seg(idx_l).cvg_dep );
                                                seg(idx_l).seq = int8( mxi-1 );
                                                seg(idx_l).ave_cvg_dep = ( ( mean( sum( seg(idx_l).cvg_dep ) )));
                                                % seg_seq_lib(idx_l,1:seg_len(idx_l)) = sub_seq_est( cvg_dep_lib(:,1:seg_len(idx_l),idx_l) );
                                                b_seg_valid(idx_s) = 0;
                                                seg_connection_mat_selected(idx_children(m), :) = 0;
                                                seg_connection_mat_selected(:, idx_children(m)) = 0;
                                            end
                                        end
                                    end
                                end
                            end
                        else
                        if (n_out_deg == 1 && n_in_deg == 0) %|| (n_in_deg == 0 && n_out_deg == 1) 
                            idx_parent = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            idx_children = find( double( ( seg_connection_mat_selected(:,idx_parent(1)))).*b_seg_valid(sidxs)' );
                            n_children = length(idx_children);
                            % check if it is one stage tree
                            b_tmp = 0;
                            for m = 1:1:n_children
                                if seg( sidxs(idx_children(m)) ).len > Min_seg_len_to_combine || b_seg_valid(sidxs(idx_children(m))) == 0
                                    b_tmp = 1;
                                end
                                if sum( double( ( seg_connection_mat_selected(idx_children(m),:))).*b_seg_valid(sidxs) ) > 1
                                    b_tmp = 1;
                                else
                                    if sum( double( ( seg_connection_mat_selected(:,idx_children(m)))).*b_seg_valid(sidxs)' ) > 0
                                        b_tmp = 1;
                                    end
                                end
                            end
                            if b_tmp == 0 && n_children > 1
                                b_tmp = 0;
                                slen_eff = zeros(n_children,1);
                                for m = 1:1:n_children
                                    idx = sidxs(idx_children(m));
                                    [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04b( sub_NumSeq2NTstr( seg(idx).seq ), 0, 0 );
                                    if tail_len_f > 0
                                        b_tmp = 1;
                                    else
                                        slen_eff(m) = seg(idx).len - tail_len_b;
                                    end
                                end
                                if b_tmp == 0
                                    [v, mx_idx] = max( slen_eff );
                                    idx_l = sidxs(idx_children(mx_idx));
                                    for m = 1:1:n_children
                                        if m ~= mx_idx
                                            idx_s = sidxs(idx_children(m));
                                            if sum( abs(seg(idx_l).seq(seg(idx_l).len-slen_eff(m)+1:seg(idx_l).len) - seg(idx_s).seq(seg(idx_s).len-slen_eff(m)+1:seg(idx_s).len)) ) > 0
                                                b_tmp = 1;
                                            end
                                        end
                                    end
                                    if b_tmp == 0
                                        for m = 1:1:n_children
                                            if m ~= mx_idx
                                                idx_s = sidxs(idx_children(m));
                                                sl = slen_eff(m); %min( seg_len(idx_l), seg_len(idx_s) );
                                                seg(idx_l).cvg_dep(:,seg(idx_l).len-sl+1:seg(idx_l).len) = ...
                                                   seg(idx_l).cvg_dep(:,seg(idx_l).len-sl+1:seg(idx_l).len) + seg(idx_s).cvg_dep(:,seg(idx_s).len-sl+1:seg(idx_s).len);
                                                seg(idx_l).seq = f03_seq_est( seg(idx_l).cvg_dep );
                                                seg(idx_l).ave_cvg_dep = ( ( mean( sum( seg(idx_l).cvg_dep ) )));
                                                b_seg_valid(idx_s) = 0;
                                                seg_connection_mat_selected(idx_children(m), :) = 0;
                                                seg_connection_mat_selected(:, idx_children(m)) = 0;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                        end
                    end
                end % for k

                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found group_idx];
                        end
                    end
                end
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            if mod( n, 20 ) == 0 || n == N_groups
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf( '/%d', n_seg_count );
            end
        end
        n_seg = n_seg_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;
        
        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
        
        if sum( b_seg_valid ) ~= length( b_seg_valid ) 
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            str_disp = sprintf(' -> %d', n_seg );
            fprintf( '%s', str_disp );
            fprintf( fp_log, '%s', str_disp );
        else
            % break;
        end

       %% remove Segs(Edges) - 2
        if kk == 1
            
        b_seg_valid = ones(1,n_seg);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        n_grp_found = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) <= 1
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    idx_r = sidxs(1);
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    if seg(idx_r).len < cfg.min_tr_length % || seg(idx_r).ave_cvg_dep < cfg.min_cvg_depth_js
                        % removed
                    else
                        n_seg_found = n_seg_found + 1;
                        n_grp_found = n_grp_found + 1;
                        seg_n(n_seg_found) = seg(sidxs(1));
                        seg_n(n_seg_found).group_index = n_grp_found;
                    end
                    n_seg_count = n_seg_count + 1;
                end
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end

                for k = 1:n_seg_selected
                    idx_r = sidxs(k); 
                    b_tmp = 0;

                    seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
                    seg(idx_r).ave_cvg_dep = ( round(mean( sum( seg(idx_r).cvg_dep )) ));

                    if seg(idx_r).ave_cvg_dep < cfg.min_cvg_depth_js
                        n_out_deg = sum( double( (seg_connection_mat_selected(k,:))) );
                        n_in_deg = sum( double( (seg_connection_mat_selected(:,k))) );
                        if n_out_deg == 0 && n_in_deg == 0
                            % b_tmp = 1;
                            if seg(idx_r).len < cfg.min_tr_length % cfg.nominal_read_length
                                % removed
                            else
                                b_tmp = 1;
                            end
                        else
                            if n_out_deg > 0 && n_in_deg > 0
                                if seg(idx_r).ave_cvg_dep < 2 && seg(idx_r).len < cfg.min_seg_length
                                    b_tmp = 1;
                                else
                                    b_tmp = 1;
                                end
                            else
                                if seg(idx_r).ave_cvg_dep < 2 && seg(idx_r).len < cfg.min_seg_length
                                    b_tmp = 1;
                                else
                                    b_tmp = 1;
                                end
                            end
                        end
                    else
                        b_tmp = 1;
                    end
                    if b_tmp == 0
                        b_seg_valid(sidxs(k)) = 0;
                    end
                end

                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                
                [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected(1:seg_cnt, 1:seg_cnt) );
                for k = 1:seg_cnt
                    seg_n(n_seg_found+k).group_index = n_grp_found + grp_idx(k);
                end
                
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found seg_n(m1+n_seg_found).group_index];
                            if seg_n(m1+n_seg_found).group_index ~= seg_n(m2+n_seg_found).group_index
                                fprintf('\n      ERROR: seg_n(m1).group_index ~= seg_n(m2).group_index (%d, %d) ', m1, m2 ); 
                            end
                        end
                    end
                end
                n_grp_found = n_grp_found + n_grps;
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            if mod( n, 20 ) == 0 || n == N_groups
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf( '/%d', n_seg_count );
            end
        end        
        n_seg = n_seg_found;
        n_grp_found_org = n_grp_found;
        N_groups = n_grp_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;

        if sum( b_seg_valid ) ~= length( b_seg_valid ) || N_groups ~= n_grp_found_org
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            % str_disp = sprintf(' >> %d', n_seg );
            str_disp = sprintf(' >>> %d(%d)', n_seg, N_groups );
            fprintf( '%s', str_disp );
            fprintf( fp_log, '%s', str_disp );
        else
        end
        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
        end
        
       %% Connect singly connected segments
        % str_disp = sprintf('   Condensing ... %d ', n_seg );
        % fprintf( '\n%s', str_disp );
        n_loop = 1;
        while(1)
            b_seg_valid = ones(1,n_seg);
            scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
            seg_n = repmat( sg, n_seg, 1 );
            seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
            seg_connection_mat_cnt_n = 0;
            n_seg_found = 0;
            n_seg_count = 0;
            for n = 1:1:N_groups
                % group_idx = group_idx_sorted(n); % Must not be sorted
                group_idx = n;
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    seg_n(1+n_seg_found) = seg(sidxs(1));
                    n_seg_found = n_seg_found + 1;
                    n_seg_count = n_seg_count + 1;
                else
                    lst_sel = find( scm_gidx == group_idx );
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    n_seg_selected = length(sidxs);
                    if isempty(lst_sel)
                        seg_connection_mat_selected = zeros(n_seg_selected);
                    else
                        tmp_idxs = zeros(n_seg,1);
                        tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                        seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                                 tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                                 ones(length(lst_sel),1), ...
                                                 n_seg_selected, n_seg_selected );
                        seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                    end
                    if n_seg_selected ~= length(sidxs)
                        fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                    end

                    for k = 1:1:n_seg_selected
                        idx_r = sidxs(k); 
                        if b_seg_valid(idx_r) > 0
                            n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            if n_out_deg == 1
                                idx_o = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                                if b_seg_valid(sidxs(idx_o)) > 0
                                    n_in_deg = sum( double( ( seg_connection_mat_selected(:,idx_o))).*b_seg_valid(sidxs)' );
                                    if n_in_deg == 1
                                        seg_len_new = double(seg(idx_r).len + seg(sidxs(idx_o)).len);
                                        seg(idx_r).cvg_dep = [seg(idx_r).cvg_dep seg(sidxs(idx_o)).cvg_dep];
                                        seg(idx_r).len = seg_len_new;
                                        [v, mxi] = max( seg(idx_r).cvg_dep );
                                        seg(idx_r).seq = int8( mxi-1 );
                                        seg(idx_r).ave_cvg_dep = ( ( mean( sum( seg(idx_r).cvg_dep ) )));

                                        seg_connection_mat_selected(k,:) = seg_connection_mat_selected(idx_o,:);
                                        seg_connection_mat_selected(idx_o,:) = 0; 

                                        b_seg_valid(sidxs(idx_o)) = 0;
                                    end
                                end
                            end
                        end
                    end
                    seg_cnt = 0;
                    for k = 1:n_seg_selected
                        if b_seg_valid(sidxs(k)) > 0
                            seg_cnt = seg_cnt + 1;
                            seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                            seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                            seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                        else
                        end
                    end
                    for m1 = 1:1:seg_cnt
                        for m2 = 1:1:seg_cnt
                            if seg_connection_mat_selected(m1,m2) > 0
                                seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                                seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                    [m1+n_seg_found m2+n_seg_found group_idx];
                            end
                        end
                    end
                    n_seg_found = n_seg_found + seg_cnt;
                    n_seg_count = n_seg_count + n_seg_selected;
                end
                if mod( n, 20 ) == 0 || n == N_groups
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                end
                Nchar = fprintf( '/%d', n_seg_count );
                end
            end
            n_seg = n_seg_found;
            seg = seg_n(1:n_seg);
            seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
            seg_connection_mat_cnt = seg_connection_mat_cnt_n;
            
            Group_Size = zeros(N_groups,1);
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
            end
            max_gsize = max( Group_Size );
            Group_Size = zeros(N_groups,1);
            Group_members = zeros( N_groups, max_gsize );
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
                Group_members(m,Group_Size(m)) = k;
            end

            if sum( b_seg_valid ) ~= length( b_seg_valid ) 
                n_loop = n_loop + 1;
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = 0;
                end
                str_disp = sprintf(' -> %d', n_seg );
                % str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
                fprintf( '%s', str_disp );
                fprintf( fp_log, '%s', str_disp );
            else
%                 if Nchar > 0
%                     fprintf(repmat('\b', 1, Nchar));
%                     Nchar = 0;
%                     fprintf(' .');
%                 end
%                 str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%                 fprintf( '%s', str_disp );
%                 fprintf( fp_log, '%s', str_disp );
                break;
            end
        end
        
    end
end
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
    Nchar = 0;
end

%% Down sizing
if b_tail_suppress > 0
    str_disp = sprintf('   Graph splitting ... %d', n_seg );
    fprintf( '\n%s', str_disp );
    % for kk = 1:1:1
    % while(1)

        b_grp_trimmed = zeros(1,n_seg);
        b_seg_valid = ones(1,n_seg);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        n_grp_found = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) <= 1
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    idx_r = sidxs(1);
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    if seg(idx_r).ave_cvg_dep < cfg.min_cvg_depth_js || seg(idx_r).len < cfg.min_tr_length
                        % removed
                    else
                        n_seg_found = n_seg_found + 1;
                        n_grp_found = n_grp_found + 1;
                        seg_n(n_seg_found) = seg(sidxs(1));
                        seg_n(n_seg_found).group_index = n_grp_found;
                    end
                    n_seg_count = n_seg_count + 1;
                end
                b_grp_trimmed(group_idx) = 1;
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end
                for k = 1:n_seg_selected
                    idx_r = sidxs(k); 
                    seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
                    seg(idx_r).ave_cvg_dep = ( round(mean( sum( seg(idx_r).cvg_dep )) ));
                end
                
%                 edges = find( seg_connection_mat_selected );
%                 num_edges = length(edges);
%                 num_eqs = sum( sum( seg_connection_mat_selected ) > 0 ) + sum( sum( seg_connection_mat_selected, 2 ) > 0 );
%                 en_mat = zeros( n_seg_selected );
%                 en_mat( edges ) = 1:1:num_edges;
%                 % num_eqs should be 2N initially
%                 if num_edges > num_eqs
%                     fprintf('\n   ERROR: num_edges > num_eqs (%d > %d) ', num_edges, num_eqs );
%                 end
%                 e_cnt = 0;
%                 y = zeros(num_eqs,1);
%                 Ts = zeros(num_eqs, num_edges);
%                 for m = 1:1:n_seg_selected
%                     % for each row, outgoing edges, corresponding segment end 
%                     if sum( seg_connection_mat_selected(m,:) ) > 0
%                         e_cnt = e_cnt + 1;
%                         y(e_cnt) = sum( seg(sidxs(m)).cvg_dep(:,end) );
%                         Ts(e_cnt, en_mat(m, en_mat(m,:) > 0 ) ) = 1;
%                     end
%                     % for each column, incoming edges, corresponding segment start
%                     if sum( seg_connection_mat_selected(:,m) ) > 0
%                         e_cnt = e_cnt + 1;
%                         y(e_cnt) = sum( seg(sidxs(m)).cvg_dep(:,1) );
%                         Ts(e_cnt, en_mat( en_mat(:,m) > 0, m ) ) = 1;
%                     end
%                 end
%                 Step_size = cfg.lasso_step_size; 
%                 N_loop = cfg.lasso_n_loops; 
%                 Lambda = 0;
%                 Wgt = ones(num_eqs,1);
%                 [edge_wgt_vec, er] = sub_tgs( y, Ts, Wgt, Lambda, Step_size, N_loop );
%                 edge_wgt_mat = zeros( n_seg_selected );
%                 edge_wgt_mat( edges ) = edge_wgt_vec;
                
                sis = 1:1:n_seg_selected;
                b_processed = zeros( n_seg_selected, 1 );
                
                while(1) 
                    [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected );
                    grp_flag = zeros( n_grps, 1 );
                    for kg = 1:1:n_grps

                        sit = sis( grp_idx == kg );
                        n_seg_sel_tmp = length(sit);
                        sc_mat_tmp = seg_connection_mat_selected( sit, sit );   
                        
                        if n_seg_sel_tmp == 1
                            grp_flag(kg) = 1;
                            b_processed( sit ) = 1;
                            % seg_connection_mat_selected( sit, sit ) = sc_mat_tmp;
                        else
                            if sum( b_processed( sit ) ) == length(sit)
                                flag = 0;
                                n_csets = 1;
                            else
                                [n_paths, path_mat, path_len, flag, nchar] = f04_path_search_v02a( sc_mat_tmp, n_max_paths, 0 );
                                Nchar = Nchar + nchar;
                                % if Nchar > 0
                                %     fprintf(repmat('\b', 1, Nchar));
                                %     Nchar = 0;
                                % end
                                n_csets = 1;
                                % if flag == 0
                                %     [n_csets, cset_size, cset_members, T, nchar] = find_cand_sets( n_paths, path_len, path_mat, ...
                                %             n_seg_sel_tmp, sc_mat_tmp, N_Csets_max, Cset_size_max, b_disp );
                                %     Nchar = Nchar + nchar;
                                % else
                                %     n_csets = -1;
                                % end
                            end
                            if flag == 0 && n_csets > 0
                                grp_flag(kg) = 1;
                                b_processed( sit ) = 1;
                            else
%                                 ew_mat_tmp = edge_wgt_mat( sit, sit );
%                                 edges_tmp = find( sc_mat_tmp );
%                                 num_edges_tmp = length(edges_tmp);
%                                 ewv_tmp = ew_mat_tmp( edges_tmp );
%                                 [ewv_tmp_sorted, ewv_idx_sorted] = sort( ewv_tmp, 'ascend' );
%                                 pcnt = 0.0001;
%                                 n_rm = ceil(pcnt*num_edges_tmp);
%                                 sc_mat_tmp( edges_tmp( ewv_idx_sorted(1:n_rm) ) ) = 0;
%                                 seg_connection_mat_selected( sit, sit ) = sc_mat_tmp;
%                                 % ew_mat_tmp( edges_tmp( ewv_idx_sorted(1:n_rm) ) ) = 0;
%                                 % edge_wgt_mat( sit, sit ) = ew_mat_tmp;

                                [axx, si] = sort( [seg(sidxs(sit)).ave_cvg_dep], 'ascend' );
                                pcnt = 0.02;
                                n_rm = ceil(pcnt*n_seg_sel_tmp);
                                b_seg_valid( sidxs( sit( si(1:n_rm) ) ) ) = 0;
                                sc_mat_tmp(si(1:n_rm),:) = 0;
                                sc_mat_tmp(:,si(1:n_rm)) = 0;
                                seg_connection_mat_selected( sit, sit ) = sc_mat_tmp;
                            end
                        end
                    end
                    if sum( grp_flag ) == n_grps
                        b_grp_trimmed(group_idx) = 1;
                        break;
                    else
                    end
                end
                
                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                
                [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected(1:seg_cnt, 1:seg_cnt) );
                for k = 1:seg_cnt
                    seg_n(n_seg_found+k).group_index = n_grp_found + grp_idx(k);
                end
                
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found seg_n(m1+n_seg_found).group_index];
                            if seg_n(m1+n_seg_found).group_index ~= seg_n(m2+n_seg_found).group_index
                                fprintf('\n      ERROR: seg_n(m1).group_index ~= seg_n(m2).group_index (%d, %d) ', m1, m2 ); 
                            end
                        end
                    end
                end
                n_grp_found = n_grp_found + n_grps;
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            % if mod( n, 20 ) == 0 || n == N_groups
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                end
                % Nchar = 0;
                Nchar = fprintf( '/%d ', n_seg_count );
            % end
        end        
        if sum( b_grp_trimmed ) == N_groups
            c_tmp = 1;
        else
            c_tmp = 0;
        end
        n_seg = n_seg_found;
        n_grp_found_org = n_grp_found;
        N_groups = n_grp_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;

        if sum( b_seg_valid ) ~= length( b_seg_valid ) || N_groups ~= n_grp_found_org
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            % str_disp = sprintf(' >>> %d', n_seg );
%             str_disp = sprintf(' >> %d(%d)', n_seg, N_groups );
%             fprintf( '%s', str_disp );
%             fprintf( fp_log, '%s', str_disp );
        else
        end
        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
               
       %% Connect singly connected segments
        % str_disp = sprintf('   Condensing ... %d ', n_seg );
        % fprintf( '\n%s', str_disp );
        n_loop = 1;
        while(1)
            b_seg_valid = ones(1,n_seg);
            scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
            seg_n = repmat( sg, n_seg, 1 );
            seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
            seg_connection_mat_cnt_n = 0;
            n_seg_found = 0;
            n_seg_count = 0;
            for n = 1:1:N_groups
                % group_idx = group_idx_sorted(n); % Must not be sorted
                group_idx = n;
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    seg_n(1+n_seg_found) = seg(sidxs(1));
                    n_seg_found = n_seg_found + 1;
                    n_seg_count = n_seg_count + 1;
                else
                    lst_sel = find( scm_gidx == group_idx );
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    n_seg_selected = length(sidxs);
                    if isempty(lst_sel)
                        seg_connection_mat_selected = zeros(n_seg_selected);
                    else
                        tmp_idxs = zeros(n_seg,1);
                        tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                        seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                                 tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                                 ones(length(lst_sel),1), ...
                                                 n_seg_selected, n_seg_selected );
                        seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                    end
                    if n_seg_selected ~= length(sidxs)
                        fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                    end

                    for k = 1:1:n_seg_selected
                        idx_r = sidxs(k); 
                        if b_seg_valid(idx_r) > 0
                            n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            if n_out_deg == 1
                                idx_o = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                                if b_seg_valid(sidxs(idx_o)) > 0
                                    n_in_deg = sum( double( ( seg_connection_mat_selected(:,idx_o))).*b_seg_valid(sidxs)' );
                                    if n_in_deg == 1
                                        seg_len_new = double(seg(idx_r).len + seg(sidxs(idx_o)).len);
                                        seg(idx_r).cvg_dep = [seg(idx_r).cvg_dep seg(sidxs(idx_o)).cvg_dep];
                                        seg(idx_r).len = seg_len_new;
                                        [v, mxi] = max( seg(idx_r).cvg_dep );
                                        seg(idx_r).seq = int8( mxi-1 );
                                        seg(idx_r).ave_cvg_dep = ( ( mean( sum( seg(idx_r).cvg_dep ) )));

                                        seg_connection_mat_selected(k,:) = seg_connection_mat_selected(idx_o,:);
                                        seg_connection_mat_selected(idx_o,:) = 0; 

                                        b_seg_valid(sidxs(idx_o)) = 0;
                                    end
                                end
                            end
                        end
                    end
                    seg_cnt = 0;
                    for k = 1:n_seg_selected
                        if b_seg_valid(sidxs(k)) > 0
                            seg_cnt = seg_cnt + 1;
                            seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                            seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                            seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                        else
                        end
                    end
                    for m1 = 1:1:seg_cnt
                        for m2 = 1:1:seg_cnt
                            if seg_connection_mat_selected(m1,m2) > 0
                                seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                                seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                    [m1+n_seg_found m2+n_seg_found group_idx];
                            end
                        end
                    end
                    n_seg_found = n_seg_found + seg_cnt;
                    n_seg_count = n_seg_count + n_seg_selected;
                end
                % if mod( n, 20 ) == 0 || n == N_groups
                    if Nchar > 0
                        fprintf(repmat('\b', 1, Nchar));
                    end
                    Nchar = fprintf( ' %d/%d ', n_seg, n_seg_count );
                % end
            end
            n_seg = n_seg_found;
            seg = seg_n(1:n_seg);
            seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
            seg_connection_mat_cnt = seg_connection_mat_cnt_n;
            
            Group_Size = zeros(N_groups,1);
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
            end
            max_gsize = max( Group_Size );
            Group_Size = zeros(N_groups,1);
            Group_members = zeros( N_groups, max_gsize );
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
                Group_members(m,Group_Size(m)) = k;
            end

            if sum( b_seg_valid ) ~= length( b_seg_valid ) 
                n_loop = n_loop + 1;
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = 0;
                end
%                 str_disp = sprintf(' -> %d', n_seg );
%                 % str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%                 fprintf( '%s', str_disp );
%                 fprintf( fp_log, '%s', str_disp );
            else
%                 if Nchar > 0
%                     fprintf(repmat('\b', 1, Nchar));
%                     Nchar = 0;
%                     fprintf(' .');
%                 end
%                 str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%                 fprintf( '%s', str_disp );
%                 fprintf( fp_log, '%s', str_disp );
                break;
            end
        end
        if c_tmp > 0
            % break;
        end
    %end
    if Nchar > 0
        fprintf(repmat('\b', 1, Nchar));
        Nchar = 0;
    end
    str_disp = sprintf(' -> %d ', n_seg );
    % str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
    fprintf( '%s', str_disp );
    fprintf( fp_log, '%s', str_disp );
end
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
    Nchar = 0;
end

%% Down sizing

if b_tail_suppress > 0
    str_disp = sprintf('   Prunning tips ..... %d', n_seg );
    fprintf( '\n%s', str_disp );
    for kk = 1:1:1

       %% remove Segs(Edges) - 3
        if kk == 1
            
        b_seg_valid = ones(1,n_seg);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        n_grp_found = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) <= 1
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    idx_r = sidxs(1);
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    if seg(idx_r).ave_cvg_dep < cfg.min_cvg_depth_js || seg(idx_r).len < cfg.min_tr_length
                        % removed
                    else
                        n_seg_found = n_seg_found + 1;
                        n_grp_found = n_grp_found + 1;
                        seg_n(n_seg_found) = seg(sidxs(1));
                        seg_n(n_seg_found).group_index = n_grp_found;
                    end
                    n_seg_count = n_seg_count + 1;
                end
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end
                for k = 1:n_seg_selected
                    idx_r = sidxs(k); 
                    seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
                    seg(idx_r).ave_cvg_dep = ( round(mean( sum( seg(idx_r).cvg_dep )) ));
                end
                
                [n_paths, path_mat, path_len, flag, nchar] = f04_path_search_v02a( seg_connection_mat_selected, n_max_paths );
                Nchar = Nchar + nchar;
                if flag == 0
                    trc_len = zeros( n_paths, 1 );
                    b_sv = zeros( n_seg_selected, 1 );
                    for k = 1:1:n_paths
                        path_idx = k;
                        if ~isempty( path_mat(path_idx,:) )
                            idx = sidxs( path_mat(path_idx,1) );
                            trc_len(k) = seg(idx).len;
                            for m = 2:1:path_len(path_idx)
                                idx = sidxs( path_mat(path_idx,m) );
                                if seg(idx).len > 0
                                    trc_len(k) = trc_len(k) + seg(idx).len;
                                end
                            end
                        end
                    end
                    mx_trc_len = max( trc_len );
                    for k = 1:1:n_paths
                        path_idx = k;
                        if ~isempty( path_mat(path_idx,:) )
                            if trc_len(path_idx) >= max( mx_trc_len*(cfg.min_rlen_percent/100), cfg.min_tr_length )
                                b_sv( path_mat(path_idx,1:path_len(path_idx)) ) = 1;
                                % for m = 1:1:path_len(path_idx)
                                %     b_sv( path_mat(path_idx,m) ) = 1;
                                % end
                            end
                        end
                    end
                    for k = 1:n_seg_selected
                        if b_sv(k) == 0 && seg(sidxs(k)).ave_cvg_dep >= 3
                            b_sv(k) = 1;
                        end
                    end
                    b_seg_valid( sidxs ) = b_sv;
                end
                
                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                
                [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected(1:seg_cnt, 1:seg_cnt) );
                for k = 1:seg_cnt
                    seg_n(n_seg_found+k).group_index = n_grp_found + grp_idx(k);
                end
                
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found seg_n(m1+n_seg_found).group_index];
                            if seg_n(m1+n_seg_found).group_index ~= seg_n(m2+n_seg_found).group_index
                                fprintf('\n      ERROR: seg_n(m1).group_index ~= seg_n(m2).group_index (%d, %d) ', m1, m2 ); 
                            end
                        end
                    end
                end
                n_grp_found = n_grp_found + n_grps;
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            if mod( n, 20 ) == 0 || n == N_groups
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf( '/%d ', n_seg_count );
            end
        end        
        n_seg = n_seg_found;
        n_grp_found_org = n_grp_found;
        N_groups = n_grp_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;

        if sum( b_seg_valid ) ~= length( b_seg_valid ) || N_groups ~= n_grp_found_org
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            str_disp = sprintf(' >>>> %d', n_seg );
            fprintf( '%s', str_disp );
            fprintf( fp_log, '%s', str_disp );
        else
        end
        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
        end
        
       %% Connect singly connected segments
        % str_disp = sprintf('   Condensing ... %d ', n_seg );
        % fprintf( '\n%s', str_disp );
        n_loop = 1;
        while(1)
            b_seg_valid = ones(1,n_seg);
            scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
            seg_n = repmat( sg, n_seg, 1 );
            seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
            seg_connection_mat_cnt_n = 0;
            n_seg_found = 0;
            n_seg_count = 0;
            for n = 1:1:N_groups
                % group_idx = group_idx_sorted(n); % Must not be sorted
                group_idx = n;
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    seg_n(1+n_seg_found) = seg(sidxs(1));
                    n_seg_found = n_seg_found + 1;
                    n_seg_count = n_seg_count + 1;
                else
                    lst_sel = find( scm_gidx == group_idx );
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    n_seg_selected = length(sidxs);
                    if isempty(lst_sel)
                        seg_connection_mat_selected = zeros(n_seg_selected);
                    else
                        tmp_idxs = zeros(n_seg,1);
                        tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                        seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                                 tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                                 ones(length(lst_sel),1), ...
                                                 n_seg_selected, n_seg_selected );
                        seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                    end
                    if n_seg_selected ~= length(sidxs)
                        fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                    end

                    for k = 1:1:n_seg_selected
                        idx_r = sidxs(k); 
                        if b_seg_valid(idx_r) > 0
                            n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            if n_out_deg == 1
                                idx_o = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                                if b_seg_valid(sidxs(idx_o)) > 0
                                    n_in_deg = sum( double( ( seg_connection_mat_selected(:,idx_o))).*b_seg_valid(sidxs)' );
                                    if n_in_deg == 1
                                        seg_len_new = double(seg(idx_r).len + seg(sidxs(idx_o)).len);
                                        seg(idx_r).cvg_dep = [seg(idx_r).cvg_dep seg(sidxs(idx_o)).cvg_dep];
                                        seg(idx_r).len = seg_len_new;
                                        [v, mxi] = max( seg(idx_r).cvg_dep );
                                        seg(idx_r).seq = int8( mxi-1 );
                                        seg(idx_r).ave_cvg_dep = ( ( mean( sum( seg(idx_r).cvg_dep ) )));

                                        seg_connection_mat_selected(k,:) = seg_connection_mat_selected(idx_o,:);
                                        seg_connection_mat_selected(idx_o,:) = 0; 

                                        b_seg_valid(sidxs(idx_o)) = 0;
                                    end
                                end
                            end
                        end
                    end
                    seg_cnt = 0;
                    for k = 1:n_seg_selected
                        if b_seg_valid(sidxs(k)) > 0
                            seg_cnt = seg_cnt + 1;
                            seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                            seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                            seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                        else
                        end
                    end
                    for m1 = 1:1:seg_cnt
                        for m2 = 1:1:seg_cnt
                            if seg_connection_mat_selected(m1,m2) > 0
                                seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                                seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                    [m1+n_seg_found m2+n_seg_found group_idx];
                            end
                        end
                    end
                    n_seg_found = n_seg_found + seg_cnt;
                    n_seg_count = n_seg_count + n_seg_selected;
                end
                if mod( n, 20 ) == 0 || n == N_groups
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                end
                Nchar = fprintf( '/%d', n_seg_count );
                end
            end
            n_seg = n_seg_found;
            seg = seg_n(1:n_seg);
            seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
            seg_connection_mat_cnt = seg_connection_mat_cnt_n;
            
            Group_Size = zeros(N_groups,1);
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
            end
            max_gsize = max( Group_Size );
            Group_Size = zeros(N_groups,1);
            Group_members = zeros( N_groups, max_gsize );
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
                Group_members(m,Group_Size(m)) = k;
            end

            if sum( b_seg_valid ) ~= length( b_seg_valid ) 
                n_loop = n_loop + 1;
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = 0;
                end
                str_disp = sprintf(' -> %d', n_seg );
                % str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
                fprintf( '%s', str_disp );
                fprintf( fp_log, '%s', str_disp );
            else
%                 if Nchar > 0
%                     fprintf(repmat('\b', 1, Nchar));
%                     Nchar = 0;
%                     fprintf(' .');
%                 end
%                 str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%                 fprintf( '%s', str_disp );
%                 fprintf( fp_log, '%s', str_disp );
                break;
            end
        end
        
    end
end
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
    Nchar = 0;
end

%% Down sizing

if b_tail_suppress > 1
    
    %% Ave. Cvg. Depth correction

    ave_cvg_dep_eff = zeros(n_seg,1);
    mlen = Edge_seg_len_min;
    Wgt_div = Cvg_dep_wgt_div; 
    e_wgt = zeros(n_seg,1);

    b_seg_valid = ones(1,n_seg);
    scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
    seg_n = repmat( sg, n_seg, 1 );
    seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
    seg_connection_mat_cnt_n = 0;
    n_seg_found = 0;
    fprintf(' ' );
    for n = 1:1:N_groups
        % group_idx = group_idx_sorted(n); % Must not be sorted
        group_idx = n;
        if Group_Size(group_idx) == 1
            sidxs = Group_members(group_idx, 1:Group_Size(group_idx));
            if length(sidxs) ~= 1
                fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
            end
            idx_r = sidxs(1);
            if seg(idx_r).len > mlen*2
                ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len-mlen) ));
                e_wgt(idx_r,1) = ceil( double(seg(idx_r).len-2*mlen)/Wgt_div );
            else
                ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len)));
                e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
            end
            n_seg_found = n_seg_found + 1;
            seg_n(n_seg_found) = seg(sidxs(1));
        else
            lst_sel = find( scm_gidx == group_idx );
            sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
            n_seg_selected = length(sidxs);
            if isempty(lst_sel)
                seg_connection_mat_selected = zeros(n_seg_selected);
            else
                tmp_idxs = zeros(n_seg,1);
                tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                         tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                         ones(length(lst_sel),1), ...
                                         n_seg_selected, n_seg_selected );
                seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
            end
            if n_seg_selected ~= length(sidxs)
                fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
            end

            for k = 1:1:n_seg_selected
                idx_r = sidxs(k);
                if seg(idx_r).len == 0
                    ave_cvg_dep_eff(idx_r) = seg(idx_r).ave_cvg_dep;
                    e_wgt(idx_r,1) = 0.2;
                    fprintf('\n      ERROR: seg(idx_r).len == 0, %d: %f ', idx_r, seg(idx_r).ave_cvg_dep );
                else
                    n_out_deg = sum( ( seg_connection_mat_selected(k,:)) );
                    n_in_deg = sum( ( seg_connection_mat_selected(:,k)) );

                    if n_out_deg == 0 && n_in_deg == 0
                        if seg(idx_r).len > mlen*2
                            ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len-mlen) ));
                            e_wgt(idx_r,1) = ceil( double(seg(idx_r).len-2*mlen)/Wgt_div );
                        else
                            ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len)));
                            e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
                        end
                        if ave_cvg_dep_eff(idx_r) < 1
                            fprintf('\n      ERROR: cdep == 0, (%d: %f) ', idx_r, seg(idx_r).ave_cvg_dep );
                        end
                    else
                        if n_out_deg == 0 || n_in_deg == 0
                            if n_out_deg == 0
                                if seg(idx_r).len > mlen
                                    ave_cvg_dep_eff(idx_r) = mean(sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len-mlen)));
                                    e_wgt(idx_r,1) = ceil( (double(seg(idx_r).len)-mlen)/Wgt_div );
                                else
                                    ave_cvg_dep_eff(idx_r) = sum( seg(idx_r).cvg_dep(:,1)); 
                                    e_wgt(idx_r,1) = 1;
                                end
                            else
                                if seg(idx_r).len > mlen
                                    ave_cvg_dep_eff(idx_r) = mean(sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len)));
                                    e_wgt(idx_r,1) = ceil( double(double(seg(idx_r).len)-mlen)/Wgt_div );
                                else
                                    ave_cvg_dep_eff(idx_r) = sum( seg(idx_r).cvg_dep(:,seg(idx_r).len)); 
                                    e_wgt(idx_r,1) = 1;
                                end
                            end
                        else
                            ave_cvg_dep_eff(idx_r) = seg(idx_r).ave_cvg_dep;
                            e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
                        end
                    end
                end
            end

            seg_cnt = 0;
            for k = 1:n_seg_selected
                if b_seg_valid(sidxs(k)) > 0
                    seg_cnt = seg_cnt + 1;
                    seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                    seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                    seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                else
                end
            end
            for m1 = 1:1:seg_cnt
                for m2 = 1:1:seg_cnt
                    if seg_connection_mat_selected(m1,m2) > 0
                        seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                        seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                            [m1+n_seg_found m2+n_seg_found group_idx];
                    end
                end
            end
            n_seg_found = n_seg_found + seg_cnt;
        end
        n_step = round(N_groups/5);
        if mod(n, n_step) == 0
            fprintf('.' );
        end
    end
    seg = seg_n(1:n_seg);
    seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
    seg_connection_mat_cnt = seg_connection_mat_cnt_n;

    %% Find groups

    n_step = round(n_seg/5);
    Group_Size = ( zeros(N_groups,1, 'uint32') );
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
        if mod(k, n_step) == 0
            fprintf('.' );
        end
    end
    Group_Size = ( zeros(N_groups,1,'uint32') );
    Group_members = ( zeros(N_groups, max(Group_Size), 'uint32') ); 
    for k = 1:1:n_seg
        m = seg(k).group_index;
        Group_Size(m) = Group_Size(m) + 1;
        Group_members( m, Group_Size(m) ) = k;
        if mod(k, n_step) == 0
            fprintf('.' );
        end
    end
    str_disp = sprintf(' N_seg: %d', n_seg );
    fprintf('%s', str_disp);
    fprintf( fp_log, '%s', str_disp );
    
    
    str_disp = sprintf('   Down sizing ....... %d', n_seg );
    fprintf( '\n%s', str_disp );
    for kk = 1:1:1

        b_seg_valid = ones(1,n_seg);
        scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
        seg_n = repmat( sg, n_seg, 1 );
        seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
        seg_connection_mat_cnt_n = 0;
        n_seg_found = 0;
        n_seg_count = 0;
        n_grp_found = 0;
        for n = 1:1:N_groups
            % group_idx = group_idx_sorted(n); % Must not be sorted
            group_idx = n;
            if Group_Size(group_idx) <= 1
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    idx_r = sidxs(1);
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    if seg(idx_r).ave_cvg_dep < cfg.min_cvg_depth_js || seg(idx_r).len < cfg.min_tr_length
                        % removed
                    else
                        n_seg_found = n_seg_found + 1;
                        n_grp_found = n_grp_found + 1;
                        seg_n(n_seg_found) = seg(sidxs(1));
                        seg_n(n_seg_found).group_index = n_grp_found;
                    end
                    n_seg_count = n_seg_count + 1;
                end
            else
                lst_sel = find( scm_gidx == group_idx );
                sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                n_seg_selected = length(sidxs);
                if isempty(lst_sel)
                    seg_connection_mat_selected = zeros(n_seg_selected);
                else
                    tmp_idxs = zeros(n_seg,1);
                    tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                    seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                             tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                             ones(length(lst_sel),1), ...
                                             n_seg_selected, n_seg_selected );
                    seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                end
                if n_seg_selected ~= length(sidxs)
                    fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                end
                for k = 1:n_seg_selected
                    idx_r = sidxs(k); 
                    seg(idx_r).seq = f03_seq_est( seg(idx_r).cvg_dep );
                    seg(idx_r).ave_cvg_dep = ( round(mean( sum( seg(idx_r).cvg_dep )) ));
                end
                
                [n_paths, path_mat, path_len, flag, nchar] = f04_path_search_v02a( seg_connection_mat_selected, n_max_paths, 0 );
                Nchar = Nchar + nchar;
                if flag == 0
                    
                    [n_csets, cset_size, cset_members, T, nchar] = find_cand_sets( n_paths, path_len, path_mat, ...
                            n_seg_selected, seg_connection_mat_selected, N_Csets_max, Cset_size_max, 1 );
                    % Nchar = Nchar + nchar;
                    if nchar > 0
                        fprintf(repmat('\b', 1, nchar));
                    end
                    
                    if n_csets < 0
                        seg_idxs = sidxs;
                        wgt = sqrt( e_wgt(seg_idxs) ); 
                        y = ave_cvg_dep_eff(seg_idxs);
                        Step_size = cfg.lasso_step_size; 
                        N_loop = cfg.lasso_n_loops;
                        Abn_thresh = 0.0001;
                        [cset_size, cset_members, ae_vec, Lambda] =  IsoLasso( y, T, wgt, Step_size, N_loop, Abn_thresh, 0 );
                        b_path_valid = zeros( n_paths, 1 );
                        b_path_valid( cset_members ) = 1;
                        b_sv = zeros( n_seg_selected, 1 );
                        
                        for k = 1:1:n_paths
                            path_idx = k;
                            if b_path_valid(path_idx) > 0 
                                for m = 1:1:path_len(path_idx)
                                    b_sv( path_mat(path_idx,m) ) = 1;
                                end
                            end
                        end
                        b_seg_valid( sidxs ) = b_sv;
                    end
                end
                
                seg_cnt = 0;
                for k = 1:n_seg_selected
                    if b_seg_valid(sidxs(k)) > 0
                        seg_cnt = seg_cnt + 1;
                        seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                        seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                        seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                    else
                    end
                end
                
                [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected(1:seg_cnt, 1:seg_cnt) );
                for k = 1:seg_cnt
                    seg_n(n_seg_found+k).group_index = n_grp_found + grp_idx(k);
                end
                
                for m1 = 1:1:seg_cnt
                    for m2 = 1:1:seg_cnt
                        if seg_connection_mat_selected(m1,m2) > 0
                            seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                            seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                [m1+n_seg_found m2+n_seg_found seg_n(m1+n_seg_found).group_index];
                            if seg_n(m1+n_seg_found).group_index ~= seg_n(m2+n_seg_found).group_index
                                fprintf('\n      ERROR: seg_n(m1).group_index ~= seg_n(m2).group_index (%d, %d) ', m1, m2 ); 
                            end
                        end
                    end
                end
                n_grp_found = n_grp_found + n_grps;
                n_seg_found = n_seg_found + seg_cnt;
                n_seg_count = n_seg_count + n_seg_selected;
            end
            % if mod( n, 20 ) == 0 || n == N_groups
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                end
                Nchar = fprintf( '/%d ', n_seg_count );
            % end
        end        
        n_seg = n_seg_found;
        n_grp_found_org = n_grp_found;
        N_groups = n_grp_found;
        seg = seg_n(1:n_seg);
        seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
        seg_connection_mat_cnt = seg_connection_mat_cnt_n;

        if sum( b_seg_valid ) ~= length( b_seg_valid ) || N_groups ~= n_grp_found_org
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
                Nchar = 0;
            end
            str_disp = sprintf(' >>>> %d', n_seg );
            fprintf( '%s', str_disp );
            fprintf( fp_log, '%s', str_disp );
        else
        end
        Group_Size = zeros(N_groups,1);
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
        end
        max_gsize = max( Group_Size );
        Group_Size = zeros(N_groups,1);
        Group_members = zeros( N_groups, max_gsize );
        for k = 1:1:n_seg
            m = seg(k).group_index;
            Group_Size(m) = Group_Size(m) + 1;
            Group_members(m,Group_Size(m)) = k;
        end
        
       %% Connect singly connected segments
        % str_disp = sprintf('   Condensing ... %d ', n_seg );
        % fprintf( '\n%s', str_disp );
        n_loop = 1;
        while(1)
            b_seg_valid = ones(1,n_seg);
            scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
            seg_n = repmat( sg, n_seg, 1 );
            seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
            seg_connection_mat_cnt_n = 0;
            n_seg_found = 0;
            n_seg_count = 0;
            for n = 1:1:N_groups
                % group_idx = group_idx_sorted(n); % Must not be sorted
                group_idx = n;
                if Group_Size(group_idx) == 1
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    if length(sidxs) ~= 1
                        fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
                    end
                    seg_n(1+n_seg_found) = seg(sidxs(1));
                    n_seg_found = n_seg_found + 1;
                    n_seg_count = n_seg_count + 1;
                else
                    lst_sel = find( scm_gidx == group_idx );
                    sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
                    n_seg_selected = length(sidxs);
                    if isempty(lst_sel)
                        seg_connection_mat_selected = zeros(n_seg_selected);
                    else
                        tmp_idxs = zeros(n_seg,1);
                        tmp_idxs( sidxs ) = (1:1:n_seg_selected);
                        seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                                 tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                                 ones(length(lst_sel),1), ...
                                                 n_seg_selected, n_seg_selected );
                        seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
                    end
                    if n_seg_selected ~= length(sidxs)
                        fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
                    end

                    for k = 1:1:n_seg_selected
                        idx_r = sidxs(k); 
                        if b_seg_valid(idx_r) > 0
                            n_out_deg = sum( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                            if n_out_deg == 1
                                idx_o = find( double( ( seg_connection_mat_selected(k,:))).*b_seg_valid(sidxs) );
                                if b_seg_valid(sidxs(idx_o)) > 0
                                    n_in_deg = sum( double( ( seg_connection_mat_selected(:,idx_o))).*b_seg_valid(sidxs)' );
                                    if n_in_deg == 1
                                        seg_len_new = double(seg(idx_r).len + seg(sidxs(idx_o)).len);
                                        seg(idx_r).cvg_dep = [seg(idx_r).cvg_dep seg(sidxs(idx_o)).cvg_dep];
                                        seg(idx_r).len = seg_len_new;
                                        [v, mxi] = max( seg(idx_r).cvg_dep );
                                        seg(idx_r).seq = int8( mxi-1 );
                                        seg(idx_r).ave_cvg_dep = ( ( mean( sum( seg(idx_r).cvg_dep ) )));

                                        seg_connection_mat_selected(k,:) = seg_connection_mat_selected(idx_o,:);
                                        seg_connection_mat_selected(idx_o,:) = 0; 

                                        b_seg_valid(sidxs(idx_o)) = 0;
                                    end
                                end
                            end
                        end
                    end
                    seg_cnt = 0;
                    for k = 1:n_seg_selected
                        if b_seg_valid(sidxs(k)) > 0
                            seg_cnt = seg_cnt + 1;
                            seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                            seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                            seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
                        else
                        end
                    end
                    for m1 = 1:1:seg_cnt
                        for m2 = 1:1:seg_cnt
                            if seg_connection_mat_selected(m1,m2) > 0
                                seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                                seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                                    [m1+n_seg_found m2+n_seg_found group_idx];
                            end
                        end
                    end
                    n_seg_found = n_seg_found + seg_cnt;
                    n_seg_count = n_seg_count + n_seg_selected;
                end
                if mod( n, 20 ) == 0 || n == N_groups
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                end
                Nchar = fprintf( '/%d', n_seg_count );
                end
            end
            n_seg = n_seg_found;
            seg = seg_n(1:n_seg);
            seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
            seg_connection_mat_cnt = seg_connection_mat_cnt_n;
            
            Group_Size = zeros(N_groups,1);
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
            end
            max_gsize = max( Group_Size );
            Group_Size = zeros(N_groups,1);
            Group_members = zeros( N_groups, max_gsize );
            for k = 1:1:n_seg
                m = seg(k).group_index;
                Group_Size(m) = Group_Size(m) + 1;
                Group_members(m,Group_Size(m)) = k;
            end

            if sum( b_seg_valid ) ~= length( b_seg_valid ) 
                n_loop = n_loop + 1;
                if Nchar > 0
                    fprintf(repmat('\b', 1, Nchar));
                    Nchar = 0;
                end
                str_disp = sprintf(' -> %d', n_seg );
                % str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
                fprintf( '%s', str_disp );
                fprintf( fp_log, '%s', str_disp );
            else
%                 if Nchar > 0
%                     fprintf(repmat('\b', 1, Nchar));
%                     Nchar = 0;
%                     fprintf(' .');
%                 end
%                 str_disp = sprintf(' >> %d(%d) ', n_seg, n_loop );
%                 fprintf( '%s', str_disp );
%                 fprintf( fp_log, '%s', str_disp );
                break;
            end
        end        
    end
end
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
    Nchar = 0;
end

%% Ave. Cvg. Depth correction

ave_cvg_dep_eff = zeros(n_seg,1);
mlen = Edge_seg_len_min;
Wgt_div = Cvg_dep_wgt_div; 
e_wgt = zeros(n_seg,1);

b_seg_valid = ones(1,n_seg);
scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';
seg_n = repmat( sg, n_seg, 1 );
seg_connection_mat_lst_n = zeros( seg_connection_mat_cnt, 3 );
seg_connection_mat_cnt_n = 0;
n_seg_found = 0;
fprintf(' ' );
for n = 1:1:N_groups
    % group_idx = group_idx_sorted(n); % Must not be sorted
    group_idx = n;
    if Group_Size(group_idx) == 1
        sidxs = Group_members(group_idx, 1:Group_Size(group_idx));
        if length(sidxs) ~= 1
            fprintf('\n      ERROR: n_seg_selected ~= 1 (%d) ', length(sidxs) );
        end
        idx_r = sidxs(1);
        if seg(idx_r).len > mlen*2
            ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len-mlen) ));
            e_wgt(idx_r,1) = ceil( double(seg(idx_r).len-2*mlen)/Wgt_div );
        else
            ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len)));
            e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
        end
        n_seg_found = n_seg_found + 1;
        seg_n(n_seg_found) = seg(sidxs(1));
    else
        lst_sel = find( scm_gidx == group_idx );
        sidxs = Group_members(group_idx, 1:Group_Size(group_idx)); 
        n_seg_selected = length(sidxs);
        if isempty(lst_sel)
            seg_connection_mat_selected = zeros(n_seg_selected);
        else
            tmp_idxs = zeros(n_seg,1);
            tmp_idxs( sidxs ) = (1:1:n_seg_selected);
            seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                     tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                     ones(length(lst_sel),1), ...
                                     n_seg_selected, n_seg_selected );
            seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
        end
        if n_seg_selected ~= length(sidxs)
            fprintf('\n      ERROR: n_seg_selected %d ~= length(sidxs) %d ', n_seg_selected, length(sidxs) );
        end

        for k = 1:1:n_seg_selected
            idx_r = sidxs(k);
            if seg(idx_r).len == 0
                ave_cvg_dep_eff(idx_r) = seg(idx_r).ave_cvg_dep;
                e_wgt(idx_r,1) = 0.2;
                fprintf('\n      ERROR: seg(idx_r).len == 0, %d: %f ', idx_r, seg(idx_r).ave_cvg_dep );
            else
                n_out_deg = sum( ( seg_connection_mat_selected(k,:)) );
                n_in_deg = sum( ( seg_connection_mat_selected(:,k)) );

                if n_out_deg == 0 && n_in_deg == 0
                    if seg(idx_r).len > mlen*2
                        ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len-mlen) ));
                        e_wgt(idx_r,1) = ceil( double(seg(idx_r).len-2*mlen)/Wgt_div );
                    else
                        ave_cvg_dep_eff(idx_r) = mean( sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len)));
                        e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
                    end
                    if ave_cvg_dep_eff(idx_r) < 1
                        fprintf('\n      ERROR: cdep == 0, (%d: %f) ', idx_r, seg(idx_r).ave_cvg_dep );
                    end
                else
                    if n_out_deg == 0 || n_in_deg == 0
                        if n_out_deg == 0
                            if seg(idx_r).len > mlen
                                ave_cvg_dep_eff(idx_r) = mean(sum( seg(idx_r).cvg_dep(:,1:seg(idx_r).len-mlen)));
                                e_wgt(idx_r,1) = ceil( (double(seg(idx_r).len)-mlen)/Wgt_div );
                            else
                                ave_cvg_dep_eff(idx_r) = sum( seg(idx_r).cvg_dep(:,1)); 
                                e_wgt(idx_r,1) = 1;
                            end
                        else
                            if seg(idx_r).len > mlen
                                ave_cvg_dep_eff(idx_r) = mean(sum( seg(idx_r).cvg_dep(:,mlen+1:seg(idx_r).len)));
                                e_wgt(idx_r,1) = ceil( double(double(seg(idx_r).len)-mlen)/Wgt_div );
                            else
                                ave_cvg_dep_eff(idx_r) = sum( seg(idx_r).cvg_dep(:,seg(idx_r).len)); 
                                e_wgt(idx_r,1) = 1;
                            end
                        end
                    else
                        ave_cvg_dep_eff(idx_r) = seg(idx_r).ave_cvg_dep;
                        e_wgt(idx_r,1) = ceil( double(seg(idx_r).len)/Wgt_div );
                    end
                end
            end
        end

        seg_cnt = 0;
        for k = 1:n_seg_selected
            if b_seg_valid(sidxs(k)) > 0
                seg_cnt = seg_cnt + 1;
                seg_n(seg_cnt+n_seg_found) = seg(sidxs(k));
                seg_connection_mat_selected(seg_cnt,:) = seg_connection_mat_selected(k,:);
                seg_connection_mat_selected(:,seg_cnt) = seg_connection_mat_selected(:,k);
            else
            end
        end
        for m1 = 1:1:seg_cnt
            for m2 = 1:1:seg_cnt
                if seg_connection_mat_selected(m1,m2) > 0
                    seg_connection_mat_cnt_n = seg_connection_mat_cnt_n + 1;
                    seg_connection_mat_lst_n( seg_connection_mat_cnt_n, : ) = ...
                        [m1+n_seg_found m2+n_seg_found group_idx];
                end
            end
        end
        n_seg_found = n_seg_found + seg_cnt;
    end
    n_step = round(N_groups/5);
    if mod(n, n_step) == 0
        fprintf('.' );
    end
end
seg = seg_n(1:n_seg);
seg_connection_mat_lst = seg_connection_mat_lst_n(1:seg_connection_mat_cnt_n,:);
seg_connection_mat_cnt = seg_connection_mat_cnt_n;

%% Find groups

n_step = round(n_seg/5);
Group_Size = ( zeros(N_groups,1, 'uint32') );
for k = 1:1:n_seg
    m = seg(k).group_index;
    Group_Size(m) = Group_Size(m) + 1;
    if mod(k, n_step) == 0
        fprintf('.' );
    end
end
Group_Size = ( zeros(N_groups,1,'uint32') );
Group_members = ( zeros(N_groups, max(Group_Size), 'uint32') ); 
for k = 1:1:n_seg
    m = seg(k).group_index;
    Group_Size(m) = Group_Size(m) + 1;
    Group_members( m, Group_Size(m) ) = k;
    if mod(k, n_step) == 0
        fprintf('.' );
    end
end
str_disp = sprintf(' N_seg: %d', n_seg );
fprintf('%s\n', str_disp);
fprintf( fp_log, '%s\n', str_disp );

%% Load reads

if R_mode > 1
    fname_rgm = sprintf('%s.rgmap', out_file_prefix );
    fp_t = fopen( fname_rgm, 'rt' );
    if fp_t < 0
        N_rgmap = 0;
        read_group_map = [];
        R_mode = 0;
    else
        N_rgmap = fscanf( fp_t, '%f', 1 );
        if N_rgmap > 0
            read_group_map = fscanf( fp_t, '%f', [3 N_rgmap] )';
        end
        fclose(fp_t);
        Nrt = max( read_group_map(:,1) );

        fprintf('   Loading reads ');
        rds_tmp.seq1 = int8(0);
        rds_tmp.seq2 = int8(0);
        rds = repmat( rds_tmp, Nrt, 1 );

        [type, fn1, ext1] = get_file_type( input_1 );
        fp_r1 = fopen( sprintf('%s.reads', fn1), 'rt' );

        [type, fn2, ext2] = get_file_type( input_2 );
        fp_r2 = fopen( sprintf('%s.reads', fn2), 'rt' );
        
        n_step = round(Nrt/20);
        for k = 1:1:Nrt
            textscan(fp_r1,'%s',1, 'Delimiter', '\n' );
            cstr = textscan(fp_r1,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr)
                rds(k).seq1 = '';
            else
                rds(k).seq1 = sub_NTstr2NumSeq( cstr{1}{:} );
            end
            textscan(fp_r2,'%s',1, 'Delimiter', '\n' );
            cstr = textscan(fp_r2,'%s',1, 'Delimiter', '\n' );
            if isempty(cstr)
                rds(k).seq2 = '';
            else
                rds(k).seq2 = sub_NTstr2NumSeq( cstr{1}{:} );
            end
            if mod(k, n_step) == 0
                fprintf('.' );
            end
        end
        fclose(fp_r1);
        fclose(fp_r2);
        fprintf(' done \n');
    end
end

%% TD and AE for each group(graph)

fname_fasta = sprintf('%s.transcriptome', fname_tr );
fp_t = fopen( fname_fasta, 'wt' );
fprintf('   Search transcript candidates ... ' );
fprintf(fp_log, '   Search transcript candidates ... ' );
scm_gidx =  [seg(seg_connection_mat_lst(1:seg_connection_mat_cnt,1)).group_index]';

Group_Size_X = zeros(N_groups,1);
for k = 1:1:n_seg
    m = seg(k).group_index;
    Group_Size_X(m) = Group_Size_X(m) + 100000 + seg(k).len*seg(k).ave_cvg_dep;
%     if mod(k, n_step) == 0
%         fprintf('.' );
%     end
end
% for k = 1:1:n_seg
%     seg(k).ave_cvg_dep = ave_cvg_dep_eff(k);
% end

[srtd, Group_idx_sorted] = sort( Group_Size_X, 'ascend' );
% Group_idx_sorted = [1:1:N_groups];

N_tr_cand = 0;
max_tr_len = 0;
Nchar = 0;
N_aligned = 0;
N_count = 0;
Group_count = 0;
N_warn = 0;
N_tr_lasso = 0;
N_tr_single = 0;
N_tr_mpc = 0;
for n = 1:1:N_groups
    
    gidx = Group_idx_sorted(n);
    n_seg_selected = double( Group_Size(gidx) );
    seg_idxs = Group_members( gidx, 1:n_seg_selected );

    if n_seg_selected <= 1
        idx = seg_idxs;
        if n_seg_selected == 1 % && seg(idx).len >= cfg.connection_threshold_short
            [tlf, tlb] = f04_check_polyA_tail_v04b( sub_NumSeq2NTstr( seg(idx).seq ), 0, 0 );
            b_tmp = 0;
            % if (seg(idx).len - max(tlf, tlb)) < cfg.connection_threshold_short || seg(idx).ave_cvg_dep < cfg.min_cvg_depth
            if seg(idx).len < cfg.min_tr_length || seg(idx).ave_cvg_dep < cfg.min_cvg_depth
                b_tmp = 1;
            end
            if b_tmp == 0
                N_tr_cand = N_tr_cand + 1;
                
                % sk = N_tr_cand;
                tpm = 1000000000*ave_cvg_dep_eff(idx)/N_bases;
                str_hdr = sprintf('>T_%d_%d_%d_%f_%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d', n, 1, Group_Size(gidx), ave_cvg_dep_eff(idx), tpm, ...
                    n, 1, seg(idx).len, (ave_cvg_dep_eff(idx)), tpm, N_reads_valid, N_bases );
                fprintf(fp_t,'%s\n', str_hdr);
                tr_seq_c = char(sub_NumSeq2NTstr( seg(idx).seq ));
                fprintf(fp_t,'%s\n', tr_seq_c);
                max_tr_len = max( max_tr_len, seg(idx).len );
                N_tr_single = N_tr_single + 1;
            end
            if ave_cvg_dep_eff(idx) < 1
                fprintf('\n      ERROR: cdep == 0, (%d: %f, %f) ', idx, seg(idx).ave_cvg_dep, ave_cvg_dep_eff(idx) );
            end
        end
        if mod( n, 20 ) == 0 || n == N_groups
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar)); 
            end
            Nchar = fprintf(' %d/%d/%d/%d ', n, gidx, N_groups, N_tr_cand );
        end
    else
        lst_sel = find( scm_gidx == gidx );
        if isempty(lst_sel)
            seg_connection_mat_selected = zeros(n_seg_selected);
        else
            tmp_idxs = zeros(n_seg,1);
            tmp_idxs( seg_idxs ) = (1:1:n_seg_selected);
            seg_connection_mat_sel = sparse( tmp_idxs(seg_connection_mat_lst(lst_sel,1)), ...
                                     tmp_idxs(seg_connection_mat_lst(lst_sel,2)), ...
                                     ones(length(lst_sel),1), ...
                                     n_seg_selected, n_seg_selected );
            seg_connection_mat_selected = sign( full( seg_connection_mat_sel ) );
        end
        if n_seg_selected ~= length(seg_idxs)
            fprintf('\n      ERROR: n_seg_selected %d ~= length(seg_idxs) %d ', n_seg_selected, length(seg_idxs) );
        end
    
        if b_disp > 0
            str_disp = sprintf('     Group %d(%d)/%d(%d):', n, gidx, N_groups, n_seg_selected );
            fprintf('\n%s', str_disp );
            fprintf(fp_log, '\n%s', str_disp );
        end
        
       %% Check loop
        % str_disp = sprintf('Checking loop ...................' );
        % fprintf('\n%s', str_disp );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Search for possible transcripts as path
        % Input:  seg_connection_mat_selected
        % Output: n_paths, path_len, path_mat
        % Nchar_save = Nchar;
        [n_paths, path_mat, path_len, flag, nchar] = f04_path_search_v02a( seg_connection_mat_selected, n_max_paths );
        if b_disp > 0
            Nchar = nchar;
        else
            Nchar = Nchar + nchar;
        end
        mn_cvg = round( min( [seg(seg_idxs).ave_cvg_dep] ) );
        if flag > 0
            % for k = 1:1:n_seg_selected
            while(1)
                [axx, si] = sort( [seg(seg_idxs).ave_cvg_dep], 'ascend' );
                bsv = ones(n_seg_selected,1);
%                 if n_seg_selected > 500
%                     pcnt = 0.4;
%                 else
%                     if n_seg_selected > 100
%                         pcnt = 0.1;
%                     else
%                         pcnt = 0.05;
%                     end
%                 end
                pcnt = 0.05;
                n_rm = ceil(pcnt*n_seg_selected);
                bsv( si(1:n_rm) ) = 0;
                mn_cvg = round(seg(seg_idxs(si(n_rm))).ave_cvg_dep);
                n_cnt = 0;
                for m = 1:n_seg_selected
                    if bsv(m) > 0
                        n_cnt = n_cnt + 1;
                        seg_idxs(n_cnt) = seg_idxs(m);
                        seg_connection_mat_selected(n_cnt,:) = seg_connection_mat_selected(m,:);
                        seg_connection_mat_selected(:,n_cnt) = seg_connection_mat_selected(:,m);
                    end
                end
                n_seg_selected = n_cnt; 
                seg_idxs = seg_idxs(1:n_seg_selected);
                seg_connection_mat_selected = seg_connection_mat_selected(1:n_seg_selected,1:n_seg_selected);
                
                [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected );
                n_p_max = 0;
                sis = [1:1:n_seg_selected];
                min_cvgs = zeros( n_grps,1 );
                flag = 0;
                for m = 1:1:n_grps
                    sit = sis( grp_idx == m );
                    sc_mat_tmp = seg_connection_mat_selected( sit, sit );
                    min_cvgs(m) = min( [ seg( seg_idxs( sit ) ).ave_cvg_dep ] );
                    [n_paths, path_mat, path_len, flag_tmp, nchar] = f04_path_search_v02a( sc_mat_tmp, n_max_paths );
                    Nchar = Nchar + nchar;
                    if flag_tmp == 0 
                        n_p_max = n_paths;
                    else
                        if flag_tmp ~= 0
                            flag = 1;
                        end
                    end
                end
                mn_cvg = max( min_cvgs );
                if flag == 0
                    break;
                else
                    if b_disp > 0
                        nchar = fprintf('%d(%d)', n_seg_selected, mn_cvg );
                        Nchar = Nchar + nchar;
                    end
                end
            end
        end
%         if (Nchar - Nchar_save) > 0 
%             fprintf(repmat('\b', 1, (Nchar - Nchar_save))); 
%         end
%         Nchar = Nchar_save;
        
        [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat_selected );
        seg_connection_mat_selected_save = seg_connection_mat_selected;
        seg_idxs_save = seg_idxs;
        n_seg_selected_save = n_seg_selected;
        sis = 1:1:n_seg_selected_save;
        for kg = 1:1:n_grps

            sit = sis( grp_idx == kg );
            n_seg_selected = length(sit);
            seg_idxs = seg_idxs_save( sit );
            sidxs = seg_idxs;
            seg_connection_mat_selected = seg_connection_mat_selected_save( sit, sit );

            if n_seg_selected == 1

                idx = seg_idxs(1);
                [tlf, tlb] = f04_check_polyA_tail_v04b( sub_NumSeq2NTstr( seg(idx).seq ), 0, 0 );
                b_tmp = 0;
                % if (seg(idx).len - max(tlf, tlb)) < cfg.connection_threshold_short || seg(idx).ave_cvg_dep < cfg.min_cvg_depth
                if seg(idx).len < cfg.min_tr_length || seg(idx).ave_cvg_dep < cfg.min_cvg_depth
                    b_tmp = 1;
                end
                if b_tmp == 0
                    N_tr_cand = N_tr_cand + 1;

                    % sk = N_tr_cand;
                    tpm = 1000000000*ave_cvg_dep_eff(idx)/N_bases;
                    str_hdr = sprintf('>T_%d_%d_%d_%f_%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d', n, 1, Group_Size(gidx), ave_cvg_dep_eff(idx), tpm, ...
                        n, 1, seg(idx).len, (ave_cvg_dep_eff(idx)), tpm, N_reads_valid, N_bases );
                    fprintf(fp_t,'%s\n', str_hdr);
                    tr_seq_c = char(sub_NumSeq2NTstr( seg(idx).seq ));
                    fprintf(fp_t,'%s\n', tr_seq_c);
                    max_tr_len = max( max_tr_len, seg(idx).len );
                    N_tr_single = N_tr_single + 1;
                end
                if ave_cvg_dep_eff(idx) < 1
                    fprintf('\n      ERROR: cdep == 0, (%d: %f, %f) ', idx, seg(idx).ave_cvg_dep, ave_cvg_dep_eff(idx) );
                end

            else

                [n_paths, path_mat, path_len, flag, nchar] = f04_path_search_v02a( seg_connection_mat_selected, n_max_paths );
                Nchar = Nchar + nchar;
                
                % Get Edge weight
                edges = find( seg_connection_mat_selected );
                num_edges = length(edges);
                num_eqs = sum( sum( seg_connection_mat_selected ) > 0 ) + sum( sum( seg_connection_mat_selected, 2 ) > 0 );
                en_mat = zeros( n_seg_selected );
                en_mat( edges ) = 1:1:num_edges;
                % num_eqs should be 2N initially
                if num_edges > num_eqs
                    fprintf('\n   ERROR: num_edges > num_eqs (%d > %d) ', num_edges, num_eqs );
                end
                e_cnt = 0;
                y = zeros(num_eqs,1);
                Ts = zeros(num_eqs, num_edges);
                for m = 1:1:n_seg_selected
                    % for each row, outgoing edges, corresponding segment end 
                    if sum( seg_connection_mat_selected(m,:) ) > 0
                        e_cnt = e_cnt + 1;
                        y(e_cnt) = sum( seg(sidxs(m)).cvg_dep(:,end) );
                        Ts(e_cnt, en_mat(m, en_mat(m,:) > 0 ) ) = 1;
                    end
                    % for each column, incoming edges, corresponding segment start
                    if sum( seg_connection_mat_selected(:,m) ) > 0
                        e_cnt = e_cnt + 1;
                        y(e_cnt) = sum( seg(sidxs(m)).cvg_dep(:,1) );
                        Ts(e_cnt, en_mat( en_mat(:,m) > 0, m ) ) = 1;
                    end
                end
                Step_size = cfg.lasso_step_size; 
                N_loop = cfg.lasso_n_loops; 
                Lambda = 0;
                Wgt = ones(num_eqs,1);
                [edge_wgt_vec, er] = sub_tgs( y, Ts, Wgt, Lambda, Step_size, N_loop );
                % edge_wgt_mat = zeros( n_seg_selected );
                % edge_wgt_mat( edges ) = edge_wgt_vec;
                
                E = zeros( num_edges, n_paths );
                for k = 1:1:n_paths
                    % for m = 1:1:path_len(k)
                    %     T( path_mat(k,m), k ) = T( path_mat(k,m), k ) + 1;
                    % end
                    for m = 2:1:path_len(k)
                        e_idx = en_mat( path_mat(k,m-1),path_mat(k,m) );
                        E( e_idx, k ) = E( e_idx, k ) + 1;
                    end
                end                

                if b_disp > 0
                    if Nchar > 0 
                        fprintf(repmat('\b', 1, Nchar)); 
                    end
                    tvec = f04_check_loop( n_seg_selected, seg_connection_mat_selected );
                    b_tmp = sum( tvec );
                    str_disp = sprintf(' N_paths ' );
                    fprintf('%s', str_disp );
                    fprintf(fp_log, '%s', str_disp );

                    str_disp = sprintf('%d(%d)', n_paths, mn_cvg );
                    if flag == 0
                        fprintf('%s', str_disp );
                        fprintf(fp_log, '%s', str_disp );
                    else
                        fprintf('%s (Max N_paths reached)', str_disp );
                        fprintf(fp_log, '%s (Max N_paths reached)', str_disp );
                    end
                    if b_tmp == 0 % if no loop detected
                        str_disp = sprintf(', ' );
                    else
                        str_disp = sprintf(', ** LOOP ** ' );
                    end
                    fprintf('%s', str_disp );
                    fprintf(fp_log, '%s', str_disp );
                end

                if n_paths > 0

                    if b_all == 0
                        % n_csets = 1;
                        cset_size(1,1) = n_paths;
                        cset_members(1,1:n_paths) = (1:1:n_paths);
                        midx = 1;
                        Abundance = zeros(n_paths,1);
                    else
                        if b_disp > 0
                            str_disp = sprintf('TD ' );
                            fprintf('%s', str_disp );
                            fprintf(fp_log, '%s', str_disp );
                        end
                        % Find CSets
                        [n_csets, cset_size, cset_members, T] = find_cand_sets( n_paths, path_len, path_mat, ...
                                n_seg_selected, seg_connection_mat_selected, N_Csets_max, Cset_size_max, b_disp );
                        
                        if b_disp > 0
                            str_disp = sprintf(' %d Csts,', n_csets );
                            fprintf('%s', str_disp );
                            fprintf(fp_log, '%s', str_disp );
                        end
                        % Output: n_csets, cset_size, cset_members
                        % n_csets: % of candidate set
                        % cset_size(1,n_csets): size(# of Transcripts) for each candidate set
                        % cset_members(n_csets,max_cset_size): node(segment) list for each candidate set
                        %    Segment index for each node should be referred to seg_idxs = Group_members(m, 1:Group_Size(m) )             

                        % cset_members(1:n_csets,1:max(cset_size(1:n_csets)))

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % select the candidate set with minimum MMSE
                        % Input: cset_size, cset_members, seg_len, ave_cvg_dep

                        if b_disp > 0
                            if n_csets ~= 1
                                str_disp = sprintf(' AE(%d) ', cset_size(1) );
                            else
                                str_disp = sprintf(' AE ' );
                            end
                            fprintf('%s', str_disp );
                            fprintf(fp_log, '%s', str_disp );
                        end

                        wgt = sqrt( e_wgt(seg_idxs) ); 
                        y = ave_cvg_dep_eff(seg_idxs);
                        Step_size = cfg.lasso_step_size;
                        N_loop = cfg.lasso_n_loops; 
                        Lambda = 0;

                        if n_csets > 1
                            % sl = double( [seg(seg_idxs).len]' );
                                                        
                            er = zeros(n_csets,1);
                            step = round(n_csets/4);
                            for k = 1:1:n_csets
                                % Build T
                                Ts = zeros( n_seg_selected, cset_size(k) );
                                Es = zeros( num_edges, cset_size(k) );
                                for l = 1:1:cset_size(k)
                                    path_idx = cset_members(k,l);
                                    Ts(:,l) = T(:,path_idx);
                                    Es(:,l) = E(:,path_idx);
                                end
                                z = [y];% edge_wgt_vec];
                                Zs = [Ts]; %Es];
                                Wgt = [wgt]; %E_cwgt.*ones(num_edges,1)];
                                [x, er_tmp] = sub_tgs( z, Zs, Wgt, Lambda, Step_size, N_loop );
                                er(k,1) = er_tmp; %(N_loop);
                                if b_disp > 0
                                    if mod(k,step) == 0
                                        fprintf( '.' );
                                    end
                                end
                            end
                            [mmse, midx] = min( er );
                            [er_srt, er_idx] = sort( er, 'ascend' );
                            b_path_valid = zeros(n_paths,1);
                            for l = 1:1:n_csets
                                if er_srt(l) <= mmse*mse_mf 
                                    csidx = er_idx(l);
                                    for m1 = 1:1:cset_size(csidx)
                                        b_path_valid( cset_members(csidx,m1) ) = 1;
                                    end
%                                     if sum(b_path_valid) > cset_size(midx) * 2;
%                                         break;
%                                     end
                                end
                            end
                            if b_disp > 0
                                fprintf(' (%d->%d)', cset_size(midx), sum(b_path_valid) );
                                fprintf(fp_log, ' (%d->%d)', cset_size(midx), sum(b_path_valid) );
                            end
                            for l = 1:1:n_paths
                                if b_path_valid(l) > 0 
                                    b_tmp = 0;
                                    for m1 = 1:1:cset_size(midx)
                                        if cset_members(midx,m1) == l
                                            b_tmp = 1;
                                        end
                                    end
                                    if b_tmp == 0
                                        cset_size(midx) = cset_size(midx) + 1;
                                        cset_members(midx,cset_size(midx)) = l;
                                    end
                                end
                            end
                            if b_disp > 0
                                if cset_size(midx) > n_seg_selected
                                    fprintf(' ** WARNING-A ** ');
                                    % N_warn = N_warn + 1;
                                end
                            end
                            N_tr_mpc = N_tr_mpc + 1;
                        else
                            if n_csets == 1
                                midx = 1;
                                N_tr_mpc = N_tr_mpc + 1;
                            else
                                z = [y]; % edge_wgt_vec];
                                Z = [T]; % E];
                                Wgt = [wgt]; % E_cwgt.*ones(num_edges,1)];
                                Step_size = cfg.lasso_step_size; 
                                N_loop = cfg.lasso_n_loops;
                                Abn_thresh = 0.0001;
                                [cset_size, cset_members, ae_vec, Lambda] =  IsoLasso( z, Z, Wgt, Step_size, N_loop, Abn_thresh, b_disp );
                                midx = 1;
                            end
                        end
                        
                        k = midx;
                        % x
                        Es = zeros( num_edges, cset_size(k) );
                        Ts = zeros( n_seg_selected, cset_size(k) );
                        for l = 1:1:cset_size(k)
                            path_idx = cset_members(k,l);
                            Ts(:,l) = T(:,path_idx);
                            Es(:,l) = E(:,path_idx);
                        end
                        z = [y]; %edge_wgt_vec];
                        Zs = [Ts]; % Es];
                        Wgt = [wgt]; % E_cwgt.*ones(num_edges,1)];
                        
                        % [Abundance, mer] = f04_Lasso_f( y, Ts, Lambda, Step_size, N_loop );
                        Step_size = cfg.lasso_step_size; 
                        N_loop = cfg.lasso_n_loops; 
                        if cset_size(midx) <= n_seg_selected || n_csets == 1
                            Lambda = 0;
                            [Abundance, er] = sub_tgs( z, Zs, Wgt, Lambda, Step_size, N_loop );
                        else
                            if n_csets < 0
                                % Lambda = 0; 
                                [Abundance, er] = sub_tgs( z, Zs, Wgt, Lambda, Step_size, N_loop );
                            else
                                Lambda = 0;
                                [Abundance, er] = sub_tgs( z, Zs, Wgt, Lambda, Step_size, N_loop );
                                N_warn = N_warn + 1;
                            end
                        end
                        if b_disp > 0
                            % fprintf(' %d', cset_size(midx) );
                            % fprintf(fp_log, ' %d', cset_size(midx) );
                        end
                    end

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Build transcript sequence for the selected set with minimum MMSE
                    % Input: midx, cset_size, cset_members, path_mat
                    %    Segment index for each node should be referred to seg_idxs = Group_members(m, 1:Group_Size(m) )
                    % str_disp = sprintf('(%d,%d). Connecting ', midx, cset_size(k) );
                    if b_disp > 0
                        % str_disp = sprintf('Connecting ');
                        % fprintf(' %s', str_disp );
                        % fprintf(fp_log, ' %s', str_disp );
                    end

                    N_tr_cand_sofar = N_tr_cand;
                    % abn_threshold = min( median( Abundance(1:cset_size(midx)) )/4, 1 );
                    trs_tmp.abn = 0;
                    trs_tmp.seq = ' ';
                    trs_tmp.len = 0;
                    trs_tmp.cvg_dep = ( zeros(4,0) );
                    trs = repmat( trs_tmp, cset_size(midx), 1 );
                    
                    for k = 1:1:cset_size(midx)
                        path_idx = cset_members(midx,k);
                        idx = seg_idxs( path_mat(path_idx,1) );
                        Tr_seq_tmp = ( seg(idx).seq );
                        for m = 2:1:path_len(path_idx)
                            idx = seg_idxs( path_mat(path_idx,m) );
                            if seg(idx).len > 0
                                seq_tmp = ( seg(idx).seq );
                                Tr_seq_tmp = [Tr_seq_tmp seq_tmp];
                            end
                        end
                        trs(k).seq = ( Tr_seq_tmp );
                        trs(k).abn = Abundance(k);
                        trs(k).len = length(trs(k).seq);
                        trs(k).cvg_dep = f03_set_cvg( Tr_seq_tmp ).*trs(k).abn;
                    end
                    
                    b_tr_valid = ones( cset_size(midx), 1 );
                    
%                     [sorted_len, si] = sort( [trs(1:end).len], 'descend' );
%                     for ks = 1:1:cset_size(midx)-1
%                         k = si(ks);
%                         if b_tr_valid(k) > 0
%                             for ms = ks+1:1:cset_size(midx)
%                                 m = si(ms);
%                                 if b_tr_valid(m) > 0 
%                                     len_l = trs(k).len;
%                                     len_s = trs(m).len;
%                                     if len_s > len_l
%                                         fprintf('\n      ERROR: len_s > len_l ');
%                                     end
%                                     Normalized_D_threshold = 0.02;
%                                     if len_s/len_l > 0.98
%                                         [pos, ref_seq_idx, dst] = sub_tco( trs(k).seq, trs(m).seq, Normalized_D_threshold, cfg.bool_ss_ind );
%                                         if pos > 0 
%                                             if pos+trs(m).len-1 > trs(k).len
%                                                 fprintf('\n      ERROR: pos+trs(m).len-1 > trs(k).len ');
%                                             end
%                                             trs(k).cvg_dep(:,pos:pos+trs(m).len-1) = trs(k).cvg_dep(:,pos:pos+trs(m).len-1) + trs(m).cvg_dep(:,1:trs(m).len);
%                                             trs(k).seq = ( f03_seq_est( trs(k).cvg_dep ) );
%                                             trs(k).abn = mean( sum( trs(k).cvg_dep ) );
%                                             b_tr_valid(m) = 0;
%                                         else
%                                             if pos < 0
%                                                 fprintf('\n      ERROR: pos < 0 ');
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
                    if b_all > 0
                        if n_csets < 0
                            N_tr_lasso = N_tr_lasso + sum( b_tr_valid );
                        end
                    end
                    
                    for k = 1:1:cset_size(midx)
                        if b_tr_valid(k) > 0
                            Tr_seq_tmp = sub_NumSeq2NTstr( trs(k).seq );
                            Abundance(k) = trs(k).abn;
                            b_tmp = 0;
                            if b_tmp == 0 && length(Tr_seq_tmp) >= cfg.connection_threshold_short
                                [tlf, tlb] = f04_check_polyA_tail_v04b( Tr_seq_tmp, 0, 0 );
                                b_tmp = 0;
                                % if (seg(idx).len - max(tlf, tlb)) < cfg.connection_threshold_short || seg(idx).ave_cvg_dep < cfg.min_cvg_depth
                                if b_all > 0 
                                    % if (length(Tr_seq_tmp) - max(tlf, tlb)) < cfg.connection_threshold_short || Abundance(k) < cfg.min_cvg_depth
                                    if length(Tr_seq_tmp) < cfg.min_tr_length || Abundance(k) < cfg.min_cvg_depth
                                        b_tmp = 1;
                                    end
                                else
                                    % if (length(Tr_seq_tmp) - max(tlf, tlb)) < cfg.connection_threshold_short
                                    if length(Tr_seq_tmp) < cfg.min_tr_length
                                        b_tmp = 1;
                                    end
                                end
                                if b_tmp == 0
                                    N_tr_cand = N_tr_cand + 1;
                                    tpm = 1000000000*Abundance(k)/N_bases;
                                    str_hdr = sprintf('>I_%d_%d_%d_%f_%f\t%d\t%d\t%d\t%f\t%f\t%d\t%d', n, k, cset_size(midx), Abundance(k), tpm, ...
                                        n, k, length(Tr_seq_tmp), Abundance(k), tpm, N_reads_valid, N_bases );
                                    fprintf(fp_t,'%s\n', str_hdr);
                                    tr_seq_c = char(Tr_seq_tmp);
                                    fprintf(fp_t,'%s\n', tr_seq_c);
                                    max_tr_len = max( max_tr_len, length(Tr_seq_tmp) );
                                end
                            end
                        end
                    end
                    N_tr_add = N_tr_cand - N_tr_cand_sofar;
                    if b_disp > 0
                        % str_disp = sprintf(' found %d trs -> %d cnds', N_tr_add, N_tr_cand );
                        str_disp = sprintf(' found %d cnds', N_tr_cand );
                        fprintf('%s', str_disp );
                        fprintf(fp_log, ' %s', str_disp );
                    end
                    % Output: N_tr_cand, Tr_len(1,N_tr_cand), Tr_seq(N_tr_cand,max_Tr_len),
                    % Tr_abn(1,N_tr_cand), Tr_gidx(1,N_tr_cand) for each group (graph)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                if b_disp == 0
                    if Nchar > 0 
                        fprintf(repmat('\b', 1, Nchar)); 
                    end
                    Nchar = fprintf(' %d/%d/%d .. Wait.', n, N_groups, N_tr_cand );
                    % Nchar = fprintf('\n %d/%d/%d/%d .. Wait. ', n, gidx, N_groups, N_tr_cand );
                end
            end
        end % kg
    end % if n_seg_selected > 1
end % for n
fclose(fp_t);
if b_disp == 0
    if Nchar > 0 
        fprintf(repmat('\b', 1, Nchar)); 
    end
    fprintf(' %d/%d/%d .. Done.', n, N_groups, N_tr_cand );
end
% fprintf(' (%5.2f)', (100.*N_aligned/N_count) );
fprintf(' (%d)', N_warn );

str_disp = sprintf('   # of Tr candidates found: %d, Max_Len: %d, N_groups: %d (%d,%d,%d) ', N_tr_cand, max_tr_len, N_groups, N_tr_single, N_tr_mpc, N_tr_lasso );
fprintf('\n%s', str_disp);
fprintf( fp_log, '\n%s', str_disp );

str_disp = sprintf('   Detected Transcripts saved to %s.', fname_fasta );
fprintf( '\n%s', str_disp );
fprintf( fp_log, '\n%s', str_disp );

fprintf('\n%s', dstr);
str_disp = sprintf(' and completed at %s', datestr(now) );
fprintf('%s\n', str_disp);

fprintf(fp_log, '\n%s', dstr );
fprintf(fp_log, '  and completed %s\n', datestr(now) );
fclose(fp_log);

end

%% %%%%%%%%% %%
%% functions %%
%% %%%%%%%%% %%

function seq_int8 = f03_seq_est( seg_cvg )
    [axx, mxi] = max( seg_cvg );
    seq_int8 = int8( mxi-1 );
end

function cvg = f03_set_cvg( seg_seq )
    seg_len = length( seg_seq );
    cvg = uint16( zeros(4,seg_len) );
    for k = 1:1:seg_len
        cvg(seg_seq(k)+1,k) = 1;
    end
end

function [n_paths, path, len, b_flag, nchar] = f04_path_search_v02a( seg_connection_mat, n_max_paths, b_disp )

    if exist('b_disp', 'var') == 0
        b_disp = 1;
    end

    n_paths = 0;
    % path = 0;
    % len = 0;
    Max_reuse_cnt = 1;
    nchar = 0;

    % find start nodes
    [nr, nc] = size(seg_connection_mat);
    b_flag = 0;
    if nr == nc
        % n_seg = nc;
        path_mat = zeros(n_max_paths, 1000);
        len_vec = zeros(n_max_paths, 1);
        b_node_start = 1 - sign( sum( seg_connection_mat ) );
        start_node_idxs = find( b_node_start );

        for k = 1:1:length(start_node_idxs)

            n_p = 1;
            n_hop = 1;
            path_tmp = zeros(n_max_paths, 1000);
            len_tmp = zeros(n_max_paths, 1);
            path_tmp(1,1) = start_node_idxs(k);
            len_tmp(1) = 1;
            % loop_chk = zeros(n_seg,1);
            n_branches = sum( seg_connection_mat(start_node_idxs(k),:) );
            if n_branches > 0
                reuse_cnt = zeros(n_max_paths,1);
                while(1)
                    b_path_valid = ones(n_max_paths,1);
                    for m = 1:1:n_p
                        if path_tmp(m,n_hop) > 0
                            n_branches = sum( seg_connection_mat(path_tmp(m,n_hop),:) );
                            if n_branches == 0
                                len_tmp(m,1) = n_hop;
                                path_tmp(m,n_hop+1) = 0;
                            else
                                idxs = find( seg_connection_mat(path_tmp(m,n_hop),:) );
                                n_p_add = -1;
                                b_loop_tmp = 0;
                                for d = 1:1:length(idxs)
                                    sm = sum(abs(sign( path_tmp(m,1:n_hop) - idxs(d) )));
                                    r_tmp = 0;
                                    if sm < n_hop
                                        r_tmp = 1;
                                    end
                                    if n_p_add < 0
                                        if reuse_cnt(m,1) + r_tmp < Max_reuse_cnt+1
                                            n_p_add = 0;
                                            reuse_cnt(m,1) = reuse_cnt(m,1) + r_tmp;
                                            path_tmp(m,n_hop+1) = idxs(d);
                                            len_tmp(m) = 0;
                                        else
                                            b_loop_tmp = b_loop_tmp + 1;
                                        end
                                    else
                                        if reuse_cnt(m,1) + r_tmp < Max_reuse_cnt+1
                                            n_p_add = n_p_add + 1;
                                            reuse_cnt(n_p+n_p_add,1) = reuse_cnt(m,1) + r_tmp;
                                            path_tmp(n_p+n_p_add,1:n_hop) = path_tmp(m,1:n_hop);
                                            path_tmp(n_p+n_p_add,n_hop+1) = idxs(d);
                                            len_tmp(n_p+n_p_add,1) = 0;
                                        else
                                            b_loop_tmp = b_loop_tmp + 1;
                                        end
                                    end
                                end
                                if n_p_add >= 0
                                    n_p = n_p + n_p_add;
                                else
                                    %if b_loop_tmp == length(idxs)
                                        b_path_valid(m) = 0;
                                    %end
                                end
                            end
                        end
                    end
                    if n_paths + n_p > n_max_paths
                        b_flag = 1;
                        break;
                    end
                    if sum( b_path_valid(1:n_p) ) < n_p
                        n_p_tmp = 0;
                        for d = 1:1:n_p
                            if b_path_valid(d) > 0
                                n_p_tmp = n_p_tmp + 1;
                                if d ~= n_p_tmp
                                    path_tmp(n_p_tmp, 1:end) = path_tmp(d, 1:end);
                                    len_tmp(n_p_tmp) = len_tmp(d);
                                    b_path_valid(n_p_tmp) = b_path_valid(d);
                                    reuse_cnt(n_p_tmp) = reuse_cnt(d);
                                end
                            end
                        end
                        n_p = n_p_tmp;
                    end
                    if sum( path_tmp(1:n_p, n_hop+1) ) == 0
                        break;
                    else
                        n_hop = n_hop + 1;
                    end
                end
            end % if n_branches > 0
            if n_paths == 0
                path_mat(1:n_p, 1:n_hop) = path_tmp(1:n_p, 1:n_hop);
                len_vec(1:n_p,1) = len_tmp(1:n_p,1);
            else
                path_mat(n_paths+1:n_paths+n_p, 1:n_hop) = path_tmp(1:n_p, 1:n_hop);
                len_vec(n_paths+1:n_paths+n_p, 1) = len_tmp(1:n_p)';
            end
            n_paths = n_paths + n_p;

            if b_disp > 0
                if mod( k, round( length(start_node_idxs)/5 ) ) == 0
                    nchar = nchar + 1;
                    fprintf( '.' );
                end
            end
        end % for k
    end
    max_path_len = max( len_vec( 1:n_paths ) );
    path = path_mat(1:n_paths,1:max_path_len);
    len = len_vec( 1:n_paths );
end

function [n_grps, grp_idx] = f04_group_search_v01a( seg_connection_mat )

    [nr, nc] = size(seg_connection_mat);
    n_grps = 0;
    grp_idx = zeros( max(nr, nc), 1 );
    
    for n = 1:1:nr
        if grp_idx( n ) == 0
            n_grps = n_grps + 1;
            grp_idx( n ) = n_grps;
            b_visit = zeros( nr, 1);
            % n_members = 1;
            while(1)
                n_new = 0;
                for k = 1:1:nr
                    if b_visit(k) == 0
                        % b_visit(k) = 1;
                        if grp_idx(k) == n_grps
                            sit1 = find( seg_connection_mat( k, : ) );
                            if ~isempty(sit1)
                            for mt = 1:1:length( sit1 )
                                m = sit1(mt);
                                if grp_idx(m) == 0 
                                    n_new = n_new + 1;
                                    grp_idx(m) = n_grps;
                                else
                                    if grp_idx(m) ~= n_grps
                                        fprintf('\n   ERROR in Group Search grp_idx(m) ~= n_grps, 1 (%d,%d) ', k, m);
                                    end
                                end
                            end
                            end
                            sit1 = find( seg_connection_mat( :, k ) );
                            if ~isempty(sit1)
                            for mt = 1:1:length( sit1 )
                                m = sit1(mt);
                                if grp_idx(m) == 0 
                                    n_new = n_new + 1;
                                    grp_idx(m) = n_grps;
                                else
                                    if grp_idx(m) ~= n_grps
                                        fprintf('\n   ERROR in Group Search grp_idx(m) ~= n_grps, 2 (%d,%d) ', k, m);
                                    end
                                end
                            end
                            end
                        end
                    end
                end
                if n_new == 0
                    break;
                else
                    % n_members = sum( grp_idx ==  n_grps );
                end
            end
        end
    end
end

%% f04_check_polyA_tail_Num_v03

function [tail_len_f, tail_len_b] = f04_check_polyA_tail_v04b( seq, Min_tail_length, d_threshold )

    Len = length(seq);
    d_th = d_threshold;
    if Len > Min_tail_length
        ts = sum( abs(sign( seq(1:Min_tail_length) - 'T' )) );
        if ts > d_th
            tail_len_b = 0;
        else
            b_tmp = 0;
            for k = 1:1:Len-Min_tail_length
                if seq(k+Min_tail_length) ~= 'T'
                    ts = ts + 1;
                end
                if ts > d_th 
                    b_tmp = 1;
                    break;
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
        if ts > d_th
            tail_len_f = 0;
        else
            for k = 1:1:Len-Min_tail_length
                if seq(Len-Min_tail_length+1-k) ~= 'A'
                    ts = ts + 1;
                end
                if ts > d_th 
                    b_tmp = 1;
                    break;
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
        if ts1 > d_th && ts2 > d_th
            tail_len_f = 0;
            tail_len_b = 0;
        else
            if ts1 <= d_th
                tail_len_f = Len;
                tail_len_b = Len;
            else
                tail_len_f = Len;
                tail_len_b = Len;
            end
        end
    end
end

function [type, fname, ext] = get_file_type( fname_ext )

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
    fprintf('Cannot open file(s) %s, ', fname_ext );
else
    d_reads = fgetl(fp);
    if d_reads(1) == '@'
        type = 2;
    else
        if d_reads(1) == '>'
            type = 1;
        else
            type = 0;
        end
    end
    fclose(fp);
end
end

function [n_csets, cset_size, cset_members, T, nchar] = find_cand_sets( n_paths, path_len, path_mat, n_seg_selected, seg_connection_mat_selected, N_Csets_max, Cset_size_max, b_disp )

    nchar = 0;
    % Find CSets
    % for kkx = 1:1:1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Build T: segment-path mapping table
        % Input:  n_paths, path_len, path_mat, n_seg_selected, 
        % Output: T, T_dst (for each group, i.e., each graph)
        Max_Weight = sum( sum( seg_connection_mat_selected ) );
        E = zeros( n_seg_selected*n_seg_selected, n_paths );
        T = zeros( n_seg_selected, n_paths );
        g_tmp = 0;
        if n_paths > 100 %n_max_paths %1000
            g_tmp = 1;
            for k = 1:1:n_paths
                for m = 1:1:path_len(k)
                    T( path_mat(k,m), k ) = T( path_mat(k,m), k ) + 1;
                end
            end
        else
            for k = 1:1:n_paths
                for m = 1:1:path_len(k)
                    T( path_mat(k,m), k ) = T( path_mat(k,m), k ) + 1;
                end
                for m = 2:1:path_len(k)
                    E( (path_mat(k,m-1)-1)*n_seg_selected+path_mat(k,m), k ) = E( (path_mat(k,m-1)-1)*n_seg_selected+path_mat(k,m), k ) + 1;
                end
            end
            T_dst = zeros( n_paths, n_paths );
            E_dst = zeros( n_paths, n_paths );
            for k = 1:1:n_paths
                for m = 1:1:n_paths
                    T_dst(k,m) = sum( abs( T(:,k)-T(:,m) ) );
                    E_dst(k,m) = sum( abs( E(:,k)-E(:,m) ) );
                end
            end
            % T_dst(k,m) indicates Hamming distance between two paths
            % i.e., # of nodes(segments) that belongs to only one path
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % search for all the candidate sets that cover the entire graph
            % Input:  T, T_dst (for each group, i.e., each graph)
            % Output: n_csets, cset_size, cset_members
            % min cover search
            N_max_csets = N_Csets_max;   % max # of csets
            N_max_size = Cset_size_max;      % max # of cset members
            n_csets = 0;           % # of csets found
            cset_size = zeros(N_max_csets,1);
            cset_members = zeros(N_max_csets,N_max_size);
            cset_members_ind = zeros(N_max_csets,n_paths);
            
            Min_n_hop = N_max_size; 
            for k = 1:1:n_paths
                % 각 path로부터 시작하여 순차 검색
                n_s = 1;    % # of csets found so far
                % n_hop = 1;  % # of hops searched so far
                cset_size_tmp = zeros(N_max_csets,1); 
                cset_members_tmp = zeros(N_max_csets,N_max_size);
                cset_members_ind_tmp = zeros(N_max_csets,n_paths);
                
                cset_size_tmp(1,1) = 1;
                cset_members_tmp(1,1) = k;
                cset_members_ind_tmp(1,k) = 1;

                x_tmp = 0;
                if n_paths == 1
                    n_csets = 1;
                    cset_size = cset_size_tmp;
                    cset_members = cset_members_tmp;
                    cset_members_ind = cset_members_ind_tmp;
                else
                    f_tmp = 0;
                    for nh = 1:1:Min_n_hop 
                        M = n_s;  % # of csets found so far
                        % 이제까지 찾은 모든 cset에 대해, 각 cset의 segment cover가
                        % 최대가 되도록 하는 path를 찾아 cset_member에 추가함. 
                        e_tmp = 0;
                        for m = 1:1:M
                            % Get segment cover, Cov, by (m)th cset
                            if cset_size_tmp(m,1) == 1 
                                Cov = abs( E(:,cset_members_tmp(m,1)) ); 
                                % cset_members_tmp(m,1): The 1st member of mth cset
                                % Cov: segments covered by cset_members_tmp(m,1)
                            else
                                %Cov = sign( sum( T(:,cset_members_tmp(m,1:cset_size_tmp(m,1))), 2 ));
                                Cov = ( sum( abs( E(:,cset_members_tmp(m,1:cset_size_tmp(m,1))) ), 2 ) > 0 );
                                % Cov: segments covered by cset_members_tmp(m,1:cset_size_tmp(m,1))
                            end
                            % Compute new cover if we select a new member
                            % for (m)th cset
                            t_dst = zeros(n_paths,1);
                            for m1 = 1:1:n_paths
                                % check if path m1 is already a cset member
                                % b_tmp = 0;
                                % for m2 = 1:1:cset_size_tmp(m,1)
                                %     if m1 == cset_members_tmp(m,m2)
                                %         b_tmp = 1;
                                %     end
                                % end
                                b_tmp = sum( m1 == cset_members_tmp(m,1:cset_size_tmp(m,1)) );
                                if b_tmp > 0 % if yes
                                    t_dst(m1,1) = 0; 
                                else % if no
                                    % The new cover if we select (m1)th path
                                    % t_dst(m1,1) = sum( sign( Cov+T(:,m1)) );
                                    t_dst(m1,1) = sum( ( Cov+abs(E(:,m1)) ) > 0 );
                                end
                            end
                            if sum( t_dst ) > 0
                                % select next member including ties
                                maxd = max( t_dst );
                                % idxs = find( 1-sign( maxd - t_dst ) ); % Cover가 최대가 되는 path의 인덱스를 찾되 tie도 포함.
                                idxs = find( ( maxd - t_dst ) == 0 ); % Cover가 최대가 되는 path의 인덱스를 찾되 tie도 포함.
                                n_branches = length(idxs); % Cover가 최대가 되는 path의 인덱스 수 
                                % cset_size_tmp와 cset_members_tmp 갱신
                                cset_size_tmp(m,1) = cset_size_tmp(m,1) + 1; 
                                cset_members_tmp(m,cset_size_tmp(m,1)) = idxs(1);
                                cset_members_ind_tmp(1,idxs(1)) = 1;
                                if n_branches > 1 
                                    cset_members_tmp(n_s+1:n_s+n_branches-1,1:cset_size_tmp(m,1)-1) = ones(n_branches-1,1)*cset_members_tmp(m,1:cset_size_tmp(m,1)-1);
                                    cset_members_tmp(n_s+1:n_s+n_branches-1,cset_size_tmp(m,1)) = idxs(2:end);
                                    cset_size_tmp(n_s+1:n_s+n_branches-1,1) = cset_size_tmp(m,1);
                                    for km = n_s+1:n_s+n_branches-1
                                        cset_members_ind_tmp(km,cset_members_tmp(km,1:cset_size_tmp(m,1))) = 1;
                                    end
                                    n_s = n_s + n_branches-1;
                                end
                            end
                            if n_csets + n_s >= N_max_csets
                                e_tmp = 1;
                                break;
                            end
                        end
                        if e_tmp > 0
                            f_tmp = 1;
                            break;
                        else
                            % check if covered: 이제까지 찾은 cset이 그래프로 연결된
                            % segment들을 모두 덮는지 확인.(즉, cover하는지 확인)
                            M = n_s;
                            b_covered = zeros(M,1);
                            for m = 1:1:M
                                %if sum( sign( sum( T(:,cset_members_tmp(m,1:cset_size_tmp(m,1))), 2 ))) == n_seg_selected
                                if sum( sign( sum( sign( E(:,cset_members_tmp(m,1:cset_size_tmp(m,1))) ), 2 ))) == Max_Weight
                                    b_covered(m,1) = 1;
                                end
                            end
                            %fprintf('%d ', n_csets);
                            % if covered
                                % Quit
                            % Otherewise (not covered)
                                % Go on
                            if sum( b_covered(1:M,1) ) >= M % > 0 
                                % 만약 하나의 cset이라도 segment들을 모두 덮으면 탐색 종료
                                Min_n_hop = nh;
                                break;
                            end
                        end
                    end % for nh = 1:1:Min_n_hop
                    if f_tmp > 0
                        g_tmp = 1;
                        break;
                    else
                        if sum( b_covered(1:n_s,1) ) > 0
                            if n_csets == 0  % 이번에 찾은 cset이 첫 cset이면
                                min_hop = nh+1 + 10;
                            else
                                min_hop = min( cset_size(1:n_csets) );
                            end

                            if min( cset_size_tmp(1:n_s,1) ) < min_hop % 이번에 찾은 cset이 첫 cset이면
                                n_csets = 0;
                                for m = 1:1:n_s
                                    if b_covered(m,1) == 1
                                        if n_csets == 0
                                            x_tmp = x_tmp + 1;
                                            n_csets = n_csets + 1;
                                            cset_size(n_csets,1) = cset_size_tmp(m,1);
                                            % cset_members_tmp_sorted = sort( cset_members_tmp(m,1:cset_size_tmp(m,1)) );
                                            % cset_members(n_csets,1:cset_size(n_csets,1)) = cset_members_tmp_sorted;
                                            cset_members(n_csets,1:cset_size(n_csets,1)) = cset_members_tmp(m,1:cset_size_tmp(m,1));
                                            cset_members_ind(n_csets, cset_members(n_csets,1:cset_size(n_csets,1))) = 1;
                                        else
                                            % cset_members_tmp_sorted = sort( cset_members_tmp(m,1:cset_size_tmp(m,1)) );
                                            t_vec = sum( abs( cset_members_ind(1:n_csets,:) - repmat( cset_members_ind_tmp(m,:), n_csets, 1 ) ), 2 );
                                            if sum(t_vec == 0) == 0
                                                x_tmp = x_tmp + 1;
                                                n_csets = n_csets + 1;
                                                cset_size(n_csets,1) = cset_size_tmp(m,1);
                                                cset_members(n_csets,1:cset_size(n_csets,1)) = cset_members_tmp(m,1:cset_size_tmp(m,1));
                                                cset_members_ind(n_csets, cset_members(n_csets,1:cset_size(n_csets,1))) = 1;
                                                % cset_members_tmp_sorted
                                            end
                                        end
                                    end
                                end
                                % str_disp = sprintf('A: k: %d, N_csets: %d', k, n_csets );
                                % disp(str_disp);
                            else % min( cset_size_tmp(1:n_s,1) ) == min_hop
                                for m = 1:1:n_s
                                    if b_covered(m,1) == 1
                                        % cset_members_tmp_sorted = sort( cset_members_tmp(m,1:cset_size_tmp(m,1)) );
                                        % check if the cset is already found
                                        t_vec = sum( abs( cset_members_ind(1:n_csets,:) - repmat( cset_members_ind_tmp(m,:), n_csets, 1 ) ), 2 );
                                        if sum(t_vec == 0) == 0
                                            x_tmp = x_tmp + 1;
                                            n_csets = n_csets + 1;
                                            cset_size(n_csets,1) = cset_size_tmp(m,1);
                                            cset_members(n_csets,1:cset_size(n_csets,1)) = cset_members_tmp(m,1:cset_size_tmp(m,1));
                                            cset_members_ind(n_csets, cset_members(n_csets,1:cset_size(n_csets,1))) = 1;
                                            % cset_members_tmp_sorted
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                %fprintf(' %d', k);
                if b_disp > 0
%                     if x_tmp > 0
%                         fprintf('*');
%                     else
                        if mod( k, round(n_paths/5) ) == 0 
                            nchar = nchar + 1;
                            fprintf('.');
                        end
%                     end
                end
            end % for k = 1:1:n_paths
        end
        if g_tmp == 0 
        else
            n_csets = -1;
            cset_size(1,1) = n_paths;
            cset_members(1,1:n_paths) = (1:1:n_paths);
        end
    % end % for kkx

end


function [cset_size, cset_members, ae_vec, Lambda] =  IsoLasso( z, Z, Wgt, Step_size, N_loop, Abn_thresh, b_disp )

    Lambda = 0.0001;
    Nlp = 40;
    N_step = 32;
    
    [x, er_tmp] = sub_tgs( z, Z, Wgt, Lambda, Step_size, N_loop );
    er_prev = er_tmp; %(N_loop);
    for mm = 1:1:Nlp 
        [x, er_tmp] = sub_tgs( z, Z, Wgt, Lambda, Step_size, N_loop );
        if er_tmp > er_prev %er_tmp(N_loop) > er_prev
            break;
        else
            er_prev = er_tmp; %(N_loop);
            Lambda = Lambda*2;
        end
        if b_disp > 0
            if mod(mm,round(Nlp/8)) == 0
                fprintf( '~' );
            end
        end
    end
    L_start = Lambda/4;
    L_stop = Lambda;
    L_step = L_stop/(N_step-1);
    er = ones(N_step,1).*10000;
    for k = 1:1:N_step
        Lambda = L_start + (k-1)*L_step;
        [x, er_tmp] = sub_tgs( z, Z, Wgt, Lambda, Step_size, N_loop );
        er(k) = er_tmp; %(N_loop);
        if k > 1
            if er(k) > er(k-1)
                break;
            end
        end
        if b_disp > 0
            if mod(k,round(N_step/8)) == 0
                fprintf( '*' );
            end
        end
    end
    [mse,sidx] = min(er);
    Lambda = L_start - (sidx-1)*L_step;
    [x, er_tmp] = sub_tgs( z, Z, Wgt, Lambda, Step_size, N_loop );
    t_vec = x >= 0; %sign( max(x-Abn_thresh,0) );
    cset_size_tmp = sum(t_vec);
    cset_members_tmp(1:cset_size_tmp) = find(t_vec);
    midx = 1;
    % n_csets = 1;
    cset_size(1) = cset_size_tmp;
    cset_members(1,1:cset_size_tmp) = cset_members_tmp(1:cset_size_tmp); % all paths are added to cset_members
    ae_vec = x;
    if b_disp > 0
        % fprintf(' %d/%d/%d/%d', mm, sidx, cset_size_tmp, n_paths );
        fprintf(' [%d->%d]', n_paths, cset_size_tmp );
        if cset_size(midx) > n_seg_selected
            fprintf(' ** WARNING-B ** ');
            % N_warn = N_warn + 1;
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
