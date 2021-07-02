clc
nb=500;
b = rand(nb,1)*45 + 5;
b=b*1e5;
N_sampled = nan(nb,1);
jj=1;
m_e = 9.10938356e-28;
c   = 2.99792458e+10;
e   = 4.803204672997660e-10; 
f1  = 0.094750;
f2  = 0.189900;
lambda2 = 1550e-8;
lambda1 = 1548e-8;
for jj=1:nb
    root=nan;
    tau = @(N) (sqrt(pi)*e^2*(10^N)*lambda1*f1/(m_e*c*b(jj)));
    g = @(N) (sqrt(1+ log((f2*lambda2)/(f1*lambda1))/log(tau(N)/log(2))));
    while isnan(root)==1
        root = fzero(@(N) (g(N)-rand-1), 14.5+ rand*2);
    end
    N_sampled(jj) = root;
end
        
figure
histogram(N_sampled)
set(get(gca, 'Title'), 'String', 'N-sampled');