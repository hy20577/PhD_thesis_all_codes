function  dydt = sing_pend(t,y)
dydt = [y(2) ; 
    -9.81*sin(y(1))];
end