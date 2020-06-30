
Cooksey_C4_detected = fitsread('data/C4_catalogs/Cooksey_C4_cat/distfiles/Cooksey_C4_detected.fits', 'binarytable');
c4_QSO_ID           = Cooksey_C4_detected{1};
c4_NCIV             = Cooksey_C4_detected{10};
c4_zCIV             = Cooksey_C4_detected{3};

f = fopen('data/C4_catalogs/Cooksey_C4_cat/processed/los_catalog','w');
for i=1:size(c4_zCIV)
    fprintf(f,'%s \n', c4_QSO_ID{i});
end

f = fopen('data/C4_catalogs/Cooksey_C4_cat/processed/c4_catalog','w');
for i=1:size(c4_zCIV)
      fprintf(f,'%s  %f %f\n', c4_QSO_ID{i}, c4_zCIV(i), c4_NCIV(i));
end

% There are some NAN valued c4_NCIV
% extract basic QSO information from Cookse_all_QSO catalog 
cooksey_all_qso_catalog = fitsread('data/dr7/distfiles/Cooksey_all_QSO.fits', 'binarytable');
all_QSO_ID              = cooksey_all_qso_catalog{1};
all_RAh                 = cooksey_all_qso_catalog{2};
all_RAm                 = cooksey_all_qso_catalog{3};
all_RAs                 = cooksey_all_qso_catalog{4};
all_DEC_Sign_d          = cooksey_all_qso_catalog{5};
all_DEC_m               = cooksey_all_qso_catalog{6};
all_DEC_s               = cooksey_all_qso_catalog{7};
all_zqso                = cooksey_all_qso_catalog{8};
all_snrs                = cooksey_all_qso_catalog{9};
all_bal_flags           = cooksey_all_qso_catalog{10};
num_quasars             = numel(all_zqso);

% Converting RA in Hour angle to ordinary degree
all_ras = all_RAh.*15 + all_RAm./4 + all_RAs./240;
all_decs = all_RAh.*15 + all_RAm./4 + all_RAs./240;

%  Converting DEC in deg-miin-sec to degree
for i=1:numel(all_DEC_Sign_d)
    if(all_DEC_Sign_d(i)>0)
        all_decs(i) = all_DEC_Sign_d(i) + all_DEC_m(i)./60 + ...
            all_DEC_s(i)./3600;
    else
        all_decs(i) = all_DEC_Sign_d(i) - all_DEC_m(i)./60 - ...
            all_DEC_s(i)./3600;
    end
end

% save catalog 
release = 'dr7';
variables_to_save = {'all_QSO_ID', 'all_ras', 'all_decs', 'all_zqso',...
    'all_snrs', 'all_bal_flags'};
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
 
% filtering bit 1: BAL
ind = (all_bal_flags==12);
filter_flags(ind) = bitset(filter_flags(ind), 2, true);
