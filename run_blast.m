
function [n_cnt, trcand, trcand_all] = run_blast( fname_ext, fname_ref, Target_cvg )

ref_q_cvg = min(Target_cvg);
[type, fname, ext] = get_file_type3( fname_ext );

fname_blst = sprintf('%s.fasta', fname );
fp = fopen( fname_blst, 'rt' );

fprintf('\nRun BLAST using Ref: %s', fname_ref );
fprintf('\n   Reading %s.fasta ...... ', fname );
Nchar = 0;
n_cnt = 0;
n_cnt2 = 0;
aline = fgets(fp);
if aline(1) == '>'
    n_cnt = 1;
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
fprintf('%d candidates ', n_cnt );

fprintf('\n   Creating Blast DB ...');
sys_command = sprintf('makeblastdb -in %s.fasta -dbtype nucl', fname );
[x, y] = system( sys_command );
% fprintf(' done ');
str_option = '-outfmt "6 qseqid sseqid qlen length slen qstart qend sstart send nident mismatch gapopen qcovs qcovhsp bitscore" -num_threads 6 -max_target_seqs 1';
fprintf(' running Blast-N ...');
sys_command = sprintf('blastn %s -db %s.fasta -query %s -out %s.tblst', str_option, fname, fname_ref, fname );
system( sys_command );
fprintf(' done ');
fn_to_delete = sprintf('%s.fasta.nhr', fname);
delete(fn_to_delete);
fn_to_delete = sprintf('%s.fasta.nin', fname);
delete(fn_to_delete);
fn_to_delete = sprintf('%s.fasta.nsq', fname);
delete(fn_to_delete);

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
trcand_org = repmat( tr_tmp, 400000, 1 );

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
        ptr = ptr + next_idx;
        [rstr, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%s', 1 );
        trcand_org( n_trcands ).sid = rstr;
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
    if trcand_org( k ).qcvg >= ref_q_cvg 
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
            end
            if b_tmp == 0
                n_trcands = n_trcands + 1;
                trcand(n_trcands) = trcand_org( k );
            end
        else
            n_trcands = n_trcands + 1;
            trcand(n_trcands) = trcand_org( k );
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

trcand = trcand(1:n_trcands);
trcand_all = trcand_all(1:n_trcands_all);

fprintf('\n   The number of Transcripts detected (Target CVG %%): ');
Cov = [trcand(1:end).qcvg];
n_trans = zeros(1, length(Target_cvg));
for k = 1:1:length(Target_cvg)
    n_trans(k) = sum( Cov >= Target_cvg(k) );
    fprintf('%d(%d), ', n_trans(k), round(Target_cvg(k)*100) );
end
fprintf(repmat('\b', 1, 2));
fprintf('\n');

% fname_out = sprintf('%s.trdtct', fname );
% fp = fopen( fname_out, 'wt' );
% for m = 1:1:n_trcands
%     fprintf( fp, '%s\t%s\t%d\t%d\t%d\t%f\t%f\t%f\t%f\n', trcand( m ).qid, trcand( m ).sid, ...
%         trcand( m ).qlen, trcand( m ).alen, trcand( m ).slen, trcand( m ).qcvg, trcand( m ).scvg, ...
%         trcand( m ).abn_est, trcand( m ).abn_true );
% end
% fclose(fp);
% fprintf('\n   Match information was saved to %s \n', fname_out );

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


