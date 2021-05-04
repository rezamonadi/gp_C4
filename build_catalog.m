% Build catalogs usable for spectra from dr16

Cooksey_C4_detected = fitsread(...
'data/C4_catalogs/Cooksey_C4_cat/distfiles/sdss_civ_cookseyetal13_update1.fit',...
'binarytable');
c4_detcted_mjd_dr7           = Cooksey_C4_detected{2};
c4_detcted_plate_dr7        = Cooksey_C4_detected{3};
c4_detcted_fiber_dr7           = Cooksey_C4_detected{4};
Z_abs_ORG             = Cooksey_C4_detected{17};
NCIV_ORG             = Cooksey_C4_detected{27};
SigmaNCIV_ORG             = Cooksey_C4_detected{28};
[nSys,dd]=size(c4_detcted_fiber_dr7);
NCIV=zeros(nSys,1);
Z_c4=zeros(nSys,1);
for i=1:nSys
    NCIV(i) = NCIV_ORG(i,1)/SigmaNCIV_ORG(i,1)^2 + NCIV_ORG(i,2)/SigmaNCIV_ORG(i,2)^2;
    NCIV(i)=NCIV(i)/(1/SigmaNCIV_ORG(i,1)^2+1/SigmaNCIV_ORG(i,2)^2);
    Z_c4(i) = (Z_abs_ORG(i,1)+Z_abs_ORG(i,2))/2;
end

c4_QSO_ID=cell(nSys,1);
f = fopen('data/C4_catalogs/Cooksey_C4_cat/processed/los_catalog','w');
for i=1:nSys

    fprintf(f,'%05i-%04i-%03i\n', c4_detcted_mjd_dr7(i), ...
    c4_detcted_plate_dr7(i), c4_detcted_fiber_dr7(i));
    c4_QSO_ID{i}=sprintf('%05i-%04i-%04i', (c4_detcted_mjd_dr7(i)), ...
    (c4_detcted_plate_dr7(i)), (c4_detcted_fiber_dr7(i)));
end



f = fopen('data/C4_catalogs/Cooksey_C4_cat/processed/c4_catalog','w');
for i=1:nSys
      fprintf(f,'%05i-%04i-%03i  %f %f\n', c4_detcted_mjd_dr7(i), ...
      c4_detcted_plate_dr7(i), c4_detcted_fiber_dr7(i), Z_c4(i), NCIV(i));
end

save('data/C4_catalogs/Cooksey_C4_cat/processed/CIV-cat.mat','c4_QSO_ID','Z_c4','NCIV');

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
cooksey_catalog = ...
fitsread('data/dr7/distfiles/dr7qso_CIV_noBAL.fit', 'binarytable');
all_plate_dr7             = cooksey_catalog{48};
all_mjd_dr7             = cooksey_catalog{47};
all_fiber_dr7             = cooksey_catalog{49};
all_RA                = cooksey_catalog{2};
all_DEC               = cooksey_catalog{3};
all_zqso                = cooksey_catalog{4};
% all_snrs                = cooksey_catalog{172};
% all_bal_flag        = cooksey_catalog{131};
num_quasars             = numel(all_zqso);

all_z_c4 = zeros(num_quasars,1);
all_z_c4 = all_z_c4 -1;
all_NCIV = zeros(num_quasars,1);
all_c4_NCIV =zeros(num_quasars,1)-1;
all_QSO_ID=cell(num_quasars,1);
for i=1:num_quasars
    all_QSO_ID{i}=sprintf('%05i-%04i-%04i', (all_mjd_dr7(i)), ...
    (all_plate_dr7(i)), (all_fiber_dr7(i)));
end
%  adding a cloumn for c4 col density if there is a c4 for a sight line
all_ind_c4 = ismember(all_QSO_ID, c4_QSO_ID);
j=0;
for i=1:num_quasars
    if all_ind_c4(i)==1
        j=j+1;
        all_z_c4(i)=Z_c4(j);
        all_c4_NCIV(i) = NCIV(j);
    end
end



% save catalog 
release = 'dr7';
variables_to_save = {'all_plate_dr7', 'all_mjd_dr7', 'all_fiber_dr7', ...
 'all_QSO_ID', 'all_RA', 'all_DEC', 'all_zqso',...
'all_z_c4', 'all_ind_c4', 'all_c4_NCIV'};
save(sprintf('%s/catalog', processed_directory(release)), ...
    variables_to_save{:}, '-v7.3');
% because QSO_IDs are cell arrays (character) we need to match them in this
% way  if(all(all_QSO_ID{1}==c4_QSO_ID{2})==true)



% to track reasons for filtering out QSOs
filter_flags = zeros(num_quasars, 1, 'uint8');
% 
% filtering bit 0: z_QSO < 1.5
ind = (all_zqso < z_qso_cut);
filter_flags(ind) = bitset(filter_flags(ind), 1, true);
 
% filtering bit 1: BAL -> the catalog is already noBAL
% filter_flags(ind) = bitset(filter_flags(:), 2, true);