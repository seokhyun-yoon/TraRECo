
function out_file_name = trareco_v063( config_file, Min_tr_abundance )

%% TraRECo main
fprintf('+-----------------------------------+\n');
fprintf('|        TraRECo version 0.63       |\n');
fprintf('+-----------------------------------+\n');
trareco_cg( config_file );
trareco_js( config_file );  
out_file_name = trareco_td( config_file );
out_file_name = trareco_tr_filter( out_file_name, Min_tr_abundance );
fprintf('+-----------------------------------+\n');
fprintf('|    Thank you for using TraRECo    |\n');
fprintf('+-----------------------------------+\n');

end


