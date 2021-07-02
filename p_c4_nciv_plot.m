filename = sprintf('%s/processed_qsos_%s', ...
    processed_directory(release), ...
    training_set_name);
load(filename);
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_zqso                = cooksey_catalog{4};
z_c13_test = all_zqso(test_ind);

z_separation = z_c13_test - map_z_c4s;





ind_has_c4 =ismember(all_QSO_ID, c4_QSO_ID) & (filter_flags==0);


DR7Q = fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL_color.fits','binarytable');
ur = DR7Q{149};
ui = DR7Q{150};
uz = DR7Q{151};
gi = DR7Q{152};
gz = DR7Q{153};
rz = DR7Q{154};


dir = sprintf('FP-dist-%s',training_set_name);
mkdir(dir);
% % odd ratio 
% odd = log_posteriors_c4 - log_posteriors_no_c4;
% x1 = odd(ind_has_c4(test_ind));
% x2 = odd(~ind_has_c4(test_ind));
% fig=figure();
% histogram(x1(x1<100),50);
% set(get(gca, 'XLabel'), 'String', 'odd(True)');
% exportgraphics(fig, sprintf('FP-dist-%s/odd-True.pdf', training_set_name), 'ContentType','vector');
% fig=figure();
% histogram(x2(x2<100),50);
% set(get(gca, 'XLabel'), 'String', 'odd(False)');
% exportgraphics(fig, sprintf('FP-dist-%s/odd-False.pdf', training_set_name), 'ContentType','vector');


% Z limitting 
x1 = p_c4(~ind_has_c4(test_ind) & map_z_c4s<2.5);
x2 = p_c4(~ind_has_c4(test_ind) & map_z_c4s>2.5 & map_z_c4s<3.5);
x3 = p_c4(~ind_has_c4(test_ind) & map_z_c4s>3.5);
fig=figure();
histogram(x1,50); 
legend('z<2.5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-z<2.5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2.5<z<3.5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2.5<z<3.5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('z>3.5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-z>3.5.pdf', training_set_name), 'ContentType','vector');


% Z-separation limitting 
x1 = p_c4(~ind_has_c4(test_ind) & z_separation<0.1);
x2 = p_c4(~ind_has_c4(test_ind) & z_separation>0.2 & z_separation<0.35);
x3 = p_c4(~ind_has_c4(test_ind) & z_separation>0.5);
fig=figure();
histogram(x1,50); 
legend('Dz<0.1')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-Dz<0.1.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('0.2<Dz<0.35')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-0.2<z<0.35.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('Dz>0.5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-z>0.5.pdf', training_set_name), 'ContentType','vector');



% Sigma limitting 
x1 = p_c4(~ind_has_c4(test_ind) & map_sigma_c4s<25e5);
x2 = p_c4(~ind_has_c4(test_ind) & map_sigma_c4s>25e5 & map_sigma_c4s<35e5);
x3 = p_c4(~ind_has_c4(test_ind) & map_sigma_c4s>35e5);
fig=figure();
histogram(x1,50); 
legend('\sigma<25km/s')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-sigma<25.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('25km/s<\sigma<35km/s')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-25<sigma<35.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('\sigma>35km/s')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-sigma>35.pdf', training_set_name), 'ContentType','vector');

%  N limitting
x1 = p_c4(~ind_has_c4(test_ind) & map_N_c4L2s<14);
x2 = p_c4(~ind_has_c4(test_ind) & map_N_c4L2s>14 & map_N_c4L2s<15);
x3 = p_c4(~ind_has_c4(test_ind) & map_N_c4L2s>15);
fig=figure();
histogram(x1,50); 
legend('N<13')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-N<13.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('13<N<14')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-13<N<14.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('N>14')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-N>14.pdf', training_set_name), 'ContentType','vector');


% u-r limitting 
x1 = p_c4(~ind_has_c4(test_ind) & ur(test_ind)<2);
x2 = p_c4(~ind_has_c4(test_ind) & ur(test_ind)>=2 & ur(test_ind)<=5);
x3 = p_c4(~ind_has_c4(test_ind) & ur(test_ind)>5);
fig=figure();
histogram(x1,50); 
legend('ur<2')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-ur<2.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2<=ur<=5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2<ur<5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('ur>5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-ur>5.pdf', training_set_name), 'ContentType','vector');


% u-i limitting 
x1 = p_c4(~ind_has_c4(test_ind) & ui(test_ind)<2);
x2 = p_c4(~ind_has_c4(test_ind) & ui(test_ind)>=2 & ui(test_ind)<=5);
x3 = p_c4(~ind_has_c4(test_ind) & ui(test_ind)>5);
fig=figure();
histogram(x1,50); 
legend('ui<2')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-ui<2.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2<ui<5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2<ui<5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('ui>5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-ui>5.pdf', training_set_name), 'ContentType','vector');

% u-z limitting 
x1 = p_c4(~ind_has_c4(test_ind) & uz(test_ind)<2);
x2 = p_c4(~ind_has_c4(test_ind) & uz(test_ind)>=2 & uz(test_ind)<=5);
x3 = p_c4(~ind_has_c4(test_ind) & uz(test_ind)>5);
fig=figure();
histogram(x1,50); 
legend('uz<2')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-uz<2.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2<=uz<=5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2<uz<5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('uz>5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-uz>5.pdf', training_set_name), 'ContentType','vector');


% g-i limitting 
x1 = p_c4(~ind_has_c4(test_ind) & gi(test_ind)<2);
x2 = p_c4(~ind_has_c4(test_ind) & gi(test_ind)>=2 & gi(test_ind)<=5);
x3 = p_c4(~ind_has_c4(test_ind) & gi(test_ind)>5);
fig=figure();
histogram(x1,50); 
legend('gi<2')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-gi<2.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2<=gi<=5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2<gi<5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('gi>5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-gi>5.pdf', training_set_name), 'ContentType','vector');



% g-z limitting 
x1 = p_c4(~ind_has_c4(test_ind) & gz(test_ind)<2);
x2 = p_c4(~ind_has_c4(test_ind) & gz(test_ind)>=2 & gz(test_ind)<=5);
x3 = p_c4(~ind_has_c4(test_ind) & gz(test_ind)>5);
fig=figure();
histogram(x1,50); 
legend('gz<2')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-gz<2.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('2<=gz<=5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-2<gz<5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('gz>5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-gz>5.pdf', training_set_name), 'ContentType','vector');




% r-z limitting 
x1 = p_c4(~ind_has_c4(test_ind) & rz(test_ind)<0.5);
x2 = p_c4(~ind_has_c4(test_ind) & rz(test_ind)>=0.5 & rz(test_ind)<=1);
x3 = p_c4(~ind_has_c4(test_ind) & rz(test_ind)>1);
fig=figure();
histogram(x1,50); 
legend('rz<0.5')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-rz<0.5.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x2,50); 
legend('0.5<=rz<=1')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-0.5<rz<1.pdf', training_set_name), 'ContentType','vector');

fig=figure();
histogram(x3,50); 
legend('rz>1')
set(get(gca, 'XLabel'), 'String', 'p-c4(False)');
exportgraphics(fig, sprintf('FP-dist-%s/p-c4(False)-rz>1.pdf', training_set_name), 'ContentType','vector');

