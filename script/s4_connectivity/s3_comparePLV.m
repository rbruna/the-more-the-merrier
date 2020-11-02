clc
clear
% close all

% Sets the paths.
config.path.corr  = '../../stats/correlations_all/';
config.path.figs  = '../../figs/violins/';
config.path.patt  = '*.mat';

% Defines the labels for the correlations to compare.
config.label1     = 'plv_1D all';
config.label2     = 'plv_99 all';
% config.label1     = 'hcoh_1D all';
% config.label2     = 'hcoh_99 all';


% Loads both sets of correlations results.
data1             = load ( sprintf ( '%s%s', config.path.corr, config.label1 ) );
data2             = load ( sprintf ( '%s%s', config.path.corr, config.label2 ) );
shape             = size ( data1.plv_inter );



% Calculates the p-value of the distribution of correlations.
[ ~, p_lab1 ]    = ttest ( data1.plv_inter ( :, : ), 0, 'tail', 'right' );
p_lab1           = reshape ( p_lab1, shape ( 2: end ) );
[ ~, p_lab2 ]    = ttest ( data2.plv_inter ( :, : ), 0, 'tail', 'right' );
p_lab2           = reshape ( p_lab2, shape ( 2: end ) );
[ ~, p_cilab1 ]  = ttest ( data1.ciplv_inter ( :, : ), 0, 'tail', 'right' );
p_cilab1         = reshape ( p_cilab1, shape ( 2: end ) );
[ ~, p_cilab2 ]  = ttest ( data2.ciplv_inter ( :, : ), 0, 'tail', 'right' );
p_cilab2         = reshape ( p_cilab2, shape ( 2: end ) );

% Calculates the difference of correlations.
[ ~, p_d ]       = ttest ( data2.plv_inter ( :, : ), data1.plv_inter ( :, : ), 'tail', 'right' );
p_d              = reshape ( p_d, shape ( 2: end ) );
[ ~, p_cid ]     = ttest ( data2.ciplv_inter ( :, : ), data1.ciplv_inter ( :, : ), 'tail', 'right' );
p_cid            = reshape ( p_cid, shape ( 2: end ) );


% % Calculates the p-values using a Wilcoxon test.
% for cindex = 1: prod ( shape ( 2: end ) )
%     
%     % Calculates the p-value of the distribution of correlations.
%     p_plv1D   ( cindex ) = signrank ( data1.plv_inter ( :, cindex ), 0, 'tail', 'right' );
%     p_plv99   ( cindex ) = signrank ( data1.plv_inter ( :, cindex ), 0, 'tail', 'right' );
%     p_ciplv1D ( cindex ) = signrank ( data1.ciplv_inter ( :, cindex ), 0, 'tail', 'right' );
%     p_ciplv99 ( cindex ) = signrank ( data1.ciplv_inter ( :, cindex ), 0, 'tail', 'right' );
%     
%     % Calculates the difference of correlations.
%     p_plv_d   ( cindex ) = signrank ( data2.plv_inter ( :, cindex ), data1.plv_inter ( :, cindex ), 'tail', 'right' );
%     p_ciplv_d ( cindex ) = signrank ( data2.ciplv_inter ( :, cindex ), data1.ciplv_inter ( :, cindex ), 'tail', 'right' );
% end




fprintf ( 1, '%s\n', config.label1 );
fprintf ( 1, '%.8f %.8f %.8f %.8f %.8f\n', min ( p_lab1 * 25, 1 ) )
fprintf ( 1, '(ci) %s\n', config.label1 );
fprintf ( 1, '%.8f %.8f %.8f %.8f %.8f\n', min ( p_cilab1 * 25, 1 ) )
fprintf ( 1, '\n' );

fprintf ( 1, '%s\n', config.label2 );
fprintf ( 1, '%.8f %.8f %.8f %.8f\n', min ( p_lab2 ( [ 1 2 3 5 ], : ) * 20, 1 ) )
fprintf ( 1, '(ci) %s\n', config.label2 );
fprintf ( 1, '%.8f %.8f %.8f %.8f\n', min ( p_cilab2 ( [ 1 2 3 5 ], : ) * 20, 1 ) )
fprintf ( 1, '\n' );

fprintf ( 1, '%s vs. %s\n', config.label1, config.label2 );
fprintf ( 1, '%.8f %.8f %.8f %.8f\n', min ( p_d ( [ 1 2 3 5 ], : ) * 20, 1 ) )
fprintf ( 1, '(ci) %s vs. %s\n', config.label1, config.label2 );
fprintf ( 1, '%.8f %.8f %.8f %.8f\n', min ( p_cid ( [ 1 2 3 5 ], : ) * 20, 1 ) )
return
hold on
plot ( data1.plv_inter ( :, 4, 2 ) )
plot ( data2.plv_inter ( :, 4, 2 ) )
