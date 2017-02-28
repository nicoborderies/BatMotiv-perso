

t0=5;
theta=1;
noise=3;

info = @(x) log(1 + (abs(x)./noise));
rt = @(x) t0 - theta.*info(x);

bounds = [-10 10];

fplot(rt,bounds);
xlabel('value');
ylabel('rating response time (s)');
