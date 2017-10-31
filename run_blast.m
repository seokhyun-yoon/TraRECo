
function [trcand, n_tr_found] = run_blast( fname_ext, fname_ref, Target_cvg, Min_Tr_Length, n_max_cands, mode )

if exist('mode', 'var') == 0
    mode = 0;
end

if exist('n_max_cands', 'var') == 0
    n_max_cands = 1;
end

if exist('Min_Tr_Length', 'var') == 0
    Min_Tr_Length = 300;
end

if exist('soption', 'var') == 0
    soption = 3;
end

ref_q_cvg = min(Target_cvg);
[type, fname, ext] = get_file_type3( fname_ext );
if isempty(ext)
    ext = 'fasta';
end
fname_tr = sprintf('%s.%s', fname, ext );
fp = fopen( fname_tr, 'rt' );

fprintf('\nRun BLAST (Changed) using Ref: %s', fname_ref );
fprintf('\n   Reading %s.fasta ...... ', fname );
Nchar = 0;
n_cnt = 0;
n_tr_found = 0;
aline = fgets(fp);
if aline(1) == '>'
    n_cnt = 0;
end
while(1)
    tline = [];
    while(1)
        bline = fgets(fp);
        if bline < 0 
            break;
        else
            if bline(1) == '>'
                aline = bline;
                break;
            else
                tline = [tline bline(1:end-1)];
            end
        end
    end
    if bline < 0 
        break;
    else
        n_cnt = n_cnt + 1;
        if length(tline) >= Min_Tr_Length
            n_tr_found = n_tr_found + 1;
        end
        if mod( n_cnt, 100 ) == 0
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf('%d', n_cnt );
        end
    end
end
fclose(fp);
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
end
fprintf('%d(%d) candidates ', n_tr_found, n_cnt );

fprintf('\n   Creating Blast DB ...');
if mod( mode, 2) == 0
    sys_command = sprintf('makeblastdb -in %s.fasta -dbtype nucl', fname );
else
    sys_command = sprintf('makeblastdb -in %s -dbtype nucl', fname_ref );
end
[x, y] = system( sys_command );
% fprintf(' done ');
str_option = sprintf( '-outfmt "6 qseqid sseqid qlen length slen qstart qend sstart send nident mismatch gapopen qcovs qcovhsp bitscore" -num_threads 4 -max_target_seqs %d', n_max_cands );
fprintf(' running Blast-N ...');
if mod( mode, 2) == 0
    sys_command = sprintf('blastn %s -db %s.fasta -query %s -out %s.tblst', str_option, fname, fname_ref, fname );
else
    sys_command = sprintf('blastn %s  -query %s.fasta -db %s -out %s.tblst', str_option, fname, fname_ref, fname );
end
system( sys_command );
fprintf(' done ');
if mod( mode, 2) == 0
    fn_to_delete = sprintf('%s.fasta.nhr', fname);
    delete(fn_to_delete);
    fn_to_delete = sprintf('%s.fasta.nin', fname);
    delete(fn_to_delete);
    fn_to_delete = sprintf('%s.fasta.nsq', fname);
    delete(fn_to_delete);
end
fname_blst = sprintf('%s.tblst', fname );
fp = fopen( fname_blst, 'rt' );

tr_tmp.qid = ' ';
tr_tmp.sid = ' ';
tr_tmp.qlen = 0;
tr_tmp.alen = 0;
tr_tmp.slen = 0;
tr_tmp.qstart = 0;
tr_tmp.qend = 0;
tr_tmp.sstart = 0;
tr_tmp.send = 0;
tr_tmp.nident = 0;
tr_tmp.nmismatch = 0;
tr_tmp.gapopen = 0;
tr_tmp.qcovs = 0;
tr_tmp.qcovhsp = 0;
tr_tmp.bitscore = 0;
tr_tmp.scovs = 0;
tr_tmp.scvg = 0;
tr_tmp.qcvg = 0;
tr_tmp.abn_true = 0;
tr_tmp.abn_est = 0;
tr_tmp.abn_rpkm_true = 0;
tr_tmp.abn_rpkm_est = 0;
tr_tmp.g_size = 0;
tr_tmp.iso_frac = 0;
tr_tmp.exp_cvg = 0;
trcand_org = repmat( tr_tmp, 1000000, 1 );

fprintf('\n   Reading tblst ... ' );
Nchar = 0;
n_trcands = 0;
n_cnt = 0;
while(1)
    aline = fgets(fp);
    if aline < 0 
        break;
    else
        n_cnt = n_cnt + 1;
        n_trcands = n_trcands + 1;
        
        ptr = 1;
        [rstr, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%s', 1 );
        trcand_org( n_trcands ).qid = rstr;
        if soption > 2
            rv = get_abn3( rstr, ':' );
            if length(rv) >= 2
                trcand_org( n_trcands ).abn_true = rv(2);
                trcand_org( n_trcands ).abn_rpkm_true = rv(1);
                if length(rv) >= 3
                    trcand_org( n_trcands ).exp_cvg = rv(3);
                end
            end
        end
        ptr = ptr + next_idx;
        [rstr, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%s', 1 );
        trcand_org( n_trcands ).sid = rstr;
        if soption > 1
            rv = get_abn3( rstr, '_' );
            if length(rv) >= 3
                trcand_org( n_trcands ).abn_est = rv(3);
                trcand_org( n_trcands ).abn_rpkm_est = rv(2);
                trcand_org( n_trcands ).g_size = rv(4); 
                trcand_org( n_trcands ).iso_frac = rv(1); 
            end
        end
        ptr = ptr + next_idx;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).qlen = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).alen = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).slen = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).qstart = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).qend = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).sstart = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).send = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).nident = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).nmismatch = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).gapopen = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).qcovs = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).qcovhsp = val;
        [val, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%e', 1 );
        ptr = ptr + next_idx;
        trcand_org( n_trcands ).bitscore = val;
        trcand_org( n_trcands ).scvg = trcand_org( n_trcands ).alen/trcand_org( n_trcands ).slen;
        trcand_org( n_trcands ).qcvg = trcand_org( n_trcands ).alen/trcand_org( n_trcands ).qlen;
        
        % tr_tmp.qid = ' ';
        % tr_tmp.sid = ' ';
        % tr_tmp.qlen = 0;
        % tr_tmp.alen = 0;
        % tr_tmp.slen = 0;
        % tr_tmp.qstart = 0;
        % tr_tmp.qend = 0;
        % tr_tmp.sstart = 0;
        % tr_tmp.send = 0;
        % tr_tmp.nident = 0;
        % tr_tmp.nmismatch = 0;
        % tr_tmp.gapopen = 0;
        % tr_tmp.qcovs = 0;
        % tr_tmp.qcovhsp = 0;
        % tr_tmp.bitscore = 0;
        % tr_tmp.scovs = 0;
        % tr_tmp.scvg = 0;
        % tr_tmp.qcvg = 0;
        % tr_tmp.abn_true = 0;
        % tr_tmp.abn_est = 0;
        % tr_tmp.abn_rpkm_true = 0;
        % tr_tmp.abn_rpkm_est = 0;
        % tr_tmp.g_size = 0;
        % tr_tmp.iso_frac = 0;
        % tr_tmp.exp_cvg = 0;
        if mod( mode, 2) == 1
            id = trcand_org( n_trcands ).qid;
            trcand_org( n_trcands ).qid = trcand_org( n_trcands ).sid;
            trcand_org( n_trcands ).sid = id;

            len = trcand_org( n_trcands ).qlen;
            trcand_org( n_trcands ).qlen = trcand_org( n_trcands ).slen;
            trcand_org( n_trcands ).slen = len;

            start = trcand_org( n_trcands ).qstart;
            trcand_org( n_trcands ).qstart = trcand_org( n_trcands ).sstart;
            trcand_org( n_trcands ).sstart = start;

            End = trcand_org( n_trcands ).qend;
            trcand_org( n_trcands ).qend = trcand_org( n_trcands ).send;
            trcand_org( n_trcands ).send = End;

            cov = trcand_org( n_trcands ).qcovs;
            trcand_org( n_trcands ).qcovs = trcand_org( n_trcands ).scovs;
            trcand_org( n_trcands ).scovs = cov;
            
            cvg = trcand_org( n_trcands ).qcvg;
            trcand_org( n_trcands ).qcvg = trcand_org( n_trcands ).scvg;
            trcand_org( n_trcands ).scvg = cvg;
        end
        
        if mod( n_cnt, 100 ) == 0
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf('%d', n_cnt );
        end
    end
end
fclose(fp);
trcand_org = trcand_org(1:n_trcands);
if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
end
fprintf('%d,', n_cnt );

Nchar = 0;
fprintf(' Selecting candidates ... ');

trcand = repmat( tr_tmp, 200000, 1 );
n_trcands = 0;
trcand_all = repmat( tr_tmp, 200000, 1 );
n_trcands_all = 0;
for k = n_cnt:-1:1
    if trcand_org( k ).qcvg >= ref_q_cvg && trcand_org( k ).slen >= Min_Tr_Length
        if n_trcands_all > 1
            b_tmp = 0;
            if b_tmp == 0
                n_trcands_all = n_trcands_all + 1;
                trcand_all(n_trcands_all) = trcand_org( k );
            end
        else
            n_trcands_all = n_trcands_all + 1;
            trcand_all(n_trcands_all) = trcand_org( k );
        end
        if mode < 2
            if n_trcands > 1
                b_tmp = 0;
                for m = 1:1:n_trcands
                    if length(trcand_org( k ).sid) == length(trcand( m ).sid) 
                        if sum( trcand_org( k ).sid == trcand( m ).sid ) == length(trcand_org( k ).sid) 
                            if trcand( m ).qcvg < trcand_org( k ).qcvg
                                trcand( m ) = trcand_org( k );
                            end
                            b_tmp = 1;
                            break;
                        end
                    end
                    if length(trcand_org( k ).qid) == length(trcand( m ).qid) 
                        if sum( trcand_org( k ).qid == trcand( m ).qid ) == length(trcand_org( k ).qid)
                            if trcand( m ).qcvg < trcand_org( k ).qcvg
                                trcand( m ) = trcand_org( k );
                            end
                            b_tmp = 1;
                            break;
                        end
                    end
                end
                if b_tmp == 0
                    n_trcands = n_trcands + 1;
                    trcand(n_trcands) = trcand_org( k );
                end
            else
                n_trcands = n_trcands + 1;
                trcand(n_trcands) = trcand_org( k );
            end
        end
    else
    end
    if mod( n_cnt-k+1, 100 ) == 0
        if Nchar > 0
            fprintf(repmat('\b', 1, Nchar));
        end
        Nchar = fprintf('%d/%d/%d', n_cnt-k+1, n_trcands_all, n_trcands );
    end
end

if mode < 2
    trcand = trcand(1:n_trcands);
else
    trcand = trcand_all(1:n_trcands_all);
    n_trcands = n_trcands_all;
end

fprintf('\n   The number of Transcripts detected (Target CVG %%): ');
Cov = [trcand(1:end).qcvg];
n_trans = zeros(1, length(Target_cvg));
for k = 1:1:length(Target_cvg)
    n_trans(k) = sum( Cov >= Target_cvg(k) );
    fprintf('%d(%d), ', n_trans(k), round(Target_cvg(k)*100) );
end
fprintf(repmat('\b', 1, 2));
fprintf('\n');

if soption > 0
    if mode == 0
        fname_out = sprintf('%s.trinfo', fname );
    else
        fname_out = sprintf('%s.trinfo%d', fname, mode );
    end
    fp = fopen( fname_out, 'wt' );
    fprintf(fp,'@ %d\n', n_tr_found);
    fprintf(fp,'> Field Name: 1. Q_ID,\t 2. S_ID,\t 3. Q_LEN,\t 4. Aligned_LEN,\t 5. S_LEN,\t 6. Q_CVG,\t 7. S_CVG,\t 8. Abndnc_Est,\t 9. Abndnc_True,\t 10. Abndnc_Est_RPKM,\t 11. Abndnc_True_RPKM,\t 12. Group_Size,\t 13. Iso_Fraction(%%)\t 14. EXPRESSED_CVG\n');
    for m = 1:1:n_trcands
        fprintf( fp, '%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%d\t%d\t%f\t%f\t%f\n', trcand( m ).qid, trcand( m ).sid, ...
            trcand( m ).qlen, trcand( m ).alen, trcand( m ).slen, trcand( m ).qcvg, trcand( m ).scvg, ...
            trcand( m ).abn_est, trcand( m ).abn_rpkm_est, trcand( m ).g_size, round(trcand( m ).iso_frac), ...
            trcand( m ).abn_true, trcand( m ).abn_rpkm_true, trcand( m ).exp_cvg );
    end
    fclose(fp);
    fprintf('Match information saved to %s \n', fname_out );
end

end

function [type, fname, ext] = get_file_type3( fname_ext )

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
end

function rv = get_abn3( t_name, c_delimiter )
    len = length(t_name);
    rv = -1.*ones(1,20);
    Ks = len;
    nrd = 0;
    for kk = 1:1:20
        for k = Ks:-1:2
            if t_name(k) == c_delimiter
                break;
            end
        end
        if k > 1
            iv = sscanf( t_name(k+1:Ks), '%f' );
            if isempty( iv )
                break;
            else
                nrd = nrd + 1;
                rv(nrd) = iv;
                Ks = k-1;
            end
        else
        end
        if k == 2
            break;
        end
    end
    if nrd > 0
        rv = rv(1:nrd);
    else
        rv = -1;
    end
end


