
cooksey_c4_catalog = fitsread('/home/reza/gp/c4/Cooksey_c4_cat/c4_catalogs/distfiles/Cooksey_C4_detected.fits', 'binarytable');
c4_QSO_ID = cooksey_c4_catalog{1};
c4_NCIV = cooksey_c4_catalog{10};
c4_zCIV = cooksey_c4_catalog{3};
f = fopen('los_catalog','w');
for i=1:size(c4_zCIV)
    
    fprintf(f,'%s \n', c4_QSO_ID{i});

end

f = fopen('c4_catalog','w');

for i=1:size(c4_zCIV)
    
    fprintf(f,'%s  %f %f\n', c4_QSO_ID{i}, c4_zCIV(i), c4_NCIV(i));

end
% There are some NAN valued c4_NCIV

% extract basic QSO information from Cookse_all_QSO catalog 
cooksey_all_qso_catalog =  fitsread('/home/reza/gp/c4/Cooksey_c4_cat/c4_catalogs/distfiles/Cooksey_all_QSO.fits', 'binarytable');
all_QSO_ID              =  cooksey_all_qso_catalog{1};
all_zqso                = cooksey_all_qso_catalog{8};
all_snrs                =  cooksey_all_qso_catalog{9};
bal_visual_flags        = cooksey_all_qso_catalog{10};
num_quasars             = numel(all_zqso);

% because QSO_IDs are cell arrays (character) we need to match them in this
% way  if(all(all_QSO_ID{1}==c4_QSO_ID{2})==true)


% % % Maybe in the future we will use them
% % % to track reasons for filtering out QSOs
% % filter_flags = zeros(num_quasars, 1, 'uint8');
% % 
% % % filtering bit 0: z_QSO < 2.15
% % ind = (z_qsos < z_qso_cut);
% % filter_flags(ind) = bitset(filter_flags(ind), 1, true);
% % 
% % % filtering bit 1: BAL
% % ind = (bal_visual_flags);
% % filter_flags(ind) = bitset(filter_flags(ind), 2, true);
