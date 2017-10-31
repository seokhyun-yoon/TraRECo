
function out_file_name = trareco_v067( config_file, Min_tr_abundance, mode )

if ~exist('mode','var')
    mode = 1;
end

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
