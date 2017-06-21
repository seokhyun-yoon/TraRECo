function trascriptome_in_fasta = trareco_tr_filter( trascriptome_in_fa, abn_th )

[type, trascriptome_in_fasta, ext] = get_file_type2( trascriptome_in_fa );

fname_blst = sprintf('%s', trascriptome_in_fa );
fp = fopen( fname_blst, 'rt' );

for k = 1:1:length(abn_th)
    if isempty( ext )
        fname_tmp = sprintf('%s_MinCvgDepth_%3.1f.fasta', trascriptome_in_fasta, round(abn_th(k)) );
    else
        fname_tmp = sprintf('%s_MinCvgDepth_%3.1f.fasta', trascriptome_in_fasta, round(abn_th(k)) );
    end
    fpwl(k) = fopen( fname_tmp, 'wt' );
end
n_rds = zeros( length(abn_th), 2 );

fprintf('Filtering .. %s .. ', trascriptome_in_fa );

disp_period = 100;
Nchar = 0;
n_cnt = 0;
while(1)
    aline = fgetl(fp);
    bline = fgetl(fp);
    if aline < 0 
        break;
    else
        if aline(1) ~= '>' || isempty(bline)
            disp(aline);
            disp(bline);
        else
        n_cnt = n_cnt + 1;
        ptr = 2;
        [name, cnt, errmsg, next_idx] = sscanf( aline(ptr:end), '%s', 1 );
        ptr = ptr + next_idx;
        val = sscanf( aline(ptr:end), '%f', 7 );
        idx = val(1);
        grp = val(2);
        abn = val(4);
        tpm = val(5);
        len = length(bline);
        nrv = val(6);
        nrb = val(7);
        seq = bline;

        for k = 1:1:length(abn_th)
            if abn < abn_th(k)
                % No action
            else
                n_rds(k,1) = n_rds(k,1) + 1;
                fprintf( fpwl(k), '>%s\t%d\t%d\t%d\t%f\t%f\t%d\t%d\n', name, idx, grp, len, abn, tpm, nrv, nrb );
                fprintf( fpwl(k), '%s\n', seq );
            end
        end
        end
        if mod( n_cnt, disp_period ) == 0
            if Nchar > 0
                fprintf(repmat('\b', 1, Nchar));
            end
            Nchar = fprintf('%d ', n_cnt );
        end
    end
end
fclose(fp);
for k = 1:1:length(abn_th)
    fclose( fpwl(k) );
end

if Nchar > 0
    fprintf(repmat('\b', 1, Nchar));
end
fprintf('%d \n', n_cnt );

fprintf('Selected Transcriptome (from %d candidates) written to \n', n_cnt );
for k = 1:1:length(abn_th)
    fname_tmp = sprintf('%s_MinCvgDepth_%3.1f.fasta', trascriptome_in_fasta, round(abn_th(k)) );
    fprintf('   Abundance >= %f -> %s  (%d transcripts) \n', ...
        abn_th(k), fname_tmp, n_rds(k) ); 
end
end

function [type, fname, ext] = get_file_type2( fname_ext )

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

