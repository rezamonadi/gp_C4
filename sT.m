clc
clear
set_parameters;
build_catalog;
generate_c4_samples;
load(sprintf('%s/civ_samples-%s', processed_directory(training_release),...
            training_set_name), 'log_nciv_samples');

m_e = 9.10938356e-28;
c   = 2.99792458e+10;
e   = 4.803204672997660e-10; 
f1  = 0.094750;
f2  = 0.189900;
lambda2 = 1550e-8;
lambda1 = 1548e-8;
for jj=1:num_C4_samples
    root=nan;
    tau = @(b) (sqrt(pi)*e^2*(10^log_nciv_samples(jj))*lambda1*f1/(m_e*c*b));
    g = @(b) (sqrt(1+ log((f2*lambda2)/(f1*lambda1))/log(tau(b)/log(2))));
    while (isnan(root)==1 | root>40 | root<5)
        root = fzero(@(b) (g(b)-rand-1), 5+35*rand);
    end

    b_sampled(jj) = root;
end
        
figure
histogram(b_sampled)
set(get(gca, 'Title'), 'String', 'b-sampled');
saveas(gcf, 'bConditionedSamplpes.pdf', 'pdf' )
figure
scatter(b_sampled, log_nciv_samples)
set(get(gca, 'XLabel'), 'String', 'b');
set(get(gca, 'YLabel'), 'String', 'N');
saveas(gcf, 'b-Conditioned-N-Samplpes.pdf', 'pdf' )
