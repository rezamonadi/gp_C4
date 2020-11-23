x = p_c4([ind_c4(test_ind) & all_c4_NCIV(test_ind)<=13]);
histogram(x)
hold on
t = p_c4([~ind_c4(test_ind) & all_c4_NCIV(test_ind)<=13]);
histogram(t)
hold off