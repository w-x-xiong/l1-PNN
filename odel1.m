function dydt = odel1(t, y, anc, c, rho, kappa, timestamps)

%t: time
%Neurons are held in y = [z', mu', lambda']' = [t0, x', d', mu', lambda']'
%anc: receiver (anchor) position matrix
%c: signal propagation speed 
%rho: augmented Lagrangian parameter
%timestamps: [t1,t2,...,tL]'

L = length(anc);

N = L + 2 + 1;

%Number of inequality constraints
K = (L^2+5*L+2)/2;

%Number of equality constraints
M = L;

dydt = zeros(N+K+M, 1);

l1approxp = 0;
for i = 1:L
    l1approxp = l1approxp + c*(exp(2*kappa*( y(3+i)+(y(1)-timestamps(i))*c) )-1)/(exp(2*kappa*(y(3+i)+(y(1)-timestamps(i))*c))+1);
end


dydt1half = -(l1approxp - y(N+1) + sum(y(N+2:N+1+L)) + c*sum(y(N+2*L+2:N+3*L+1)) + 2*c*sum(y(N+3*L+2:N+K)) + rho*(y(N+1)^2*y(1) + sum((y(N+2:N+1+L).^2).*(ones(L,1)*y(1)-timestamps)) + c*sum((y(N+2*L+2:N+3*L+1).^2).*(y(3+1:3+L) - (timestamps-ones(L,1)*y(1))*c)) ));
lsttm = 0;
for i = 1:L-1
   for j = i+1:L
       lsttm = lsttm + 2*c*y(N + (2*L-i)*(i-1)/2 + j - i + 3*L + 1)^2*((2*y(1)-timestamps(i)-timestamps(j))*c + norm(anc(:,i) - anc(:,j)));
   end
end
%dto/dt
dydt(1) = dydt1half - rho*lsttm;



%dx/dt
for i = 1:L
dydt(2:3) = dydt(2:3) - 2*(y(N+K+i) + rho*y(N+K+i)^2*(y(3+i)^2 - norm(y(2:3) - anc(:,i))^2))*(anc(:,i) - y(2:3));
end


%ddi/dt
for i = 1:L
    dydt(3+i) = -( (exp(2*kappa*(y(3+i)+(y(1)-timestamps(i))*c))-1)/(exp(2*kappa*(y(3+i)+(y(1)-timestamps(i))*c))+1) - y(N+L+1+i) + y(N+i+2*L+1) + 2*y(N+K+i)*y(3+i) + rho*( y(N+L+1+i)^2*y(3+i) + y(N+i+2*L+1)^2*(y(3+i) - (timestamps(i)-y(1))*c) + 2*y(N+K+i)^2*y(3+i)*(y(3+i)^2 - norm(y(2:3) - anc(:,i))^2) ) );
end



dydt(N+1) = -y(N+1) + max(0, y(N+1) - y(1));

for i = 1:L
dydt(N+1+i) = -y(N+1+i) + max(0, y(N+1+i) + y(1) - timestamps(i));
dydt(N+1+i+L) = -y(N+1+i+L) + max(0, y(N+1+i+L)-y(3+i));
dydt(N+1+i+2*L) = -y(N+1+i+2*L) + max(0, y(N+1+i+2*L)+y(3+i) - (timestamps(i) - y(1))*c);
end

for i = 1:L-1
    for j = i+1:L
        dydt(N+1+3*L+(2*L-i)*(i-1)/2 + j - i) = -y(N+1+3*L+(2*L-i)*(i-1)/2 + j - i) + max(0, y(N+1+3*L+(2*L-i)*(i-1)/2 + j - i)+ ((2*y(1)-timestamps(i)-timestamps(j))*c + norm(anc(:,i) - anc(:,j))) );
    end
end

for i = 1:L
    dydt(N+K+i) = (y(3+i)^2 - norm(y(2:3) - anc(:,i))^2);
end

end

