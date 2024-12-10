function Freq = deterministicKuramoto_freq(phases, A, K)
y= phases;
M = size(phases,2);
n = size(phases, 1);
Freq = zeros(size(phases));
w = rand(1,M)*2*pi-pi;   %internal frequencies has 0 mean, [-pi, pi].

for k=1:n
    for i=1:M
        sum1 = 0;
        sum2 = 0;
        for j=1:M
            sum1 = sum1+ A(i,j)*sin(y(k,j)-y(k,i));
        end
    theta(i) = w(i) + (K/M)*sum1;  
    end
Freq(k,:) = theta;
end

end

