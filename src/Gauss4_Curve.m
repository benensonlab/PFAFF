function [y] = Gauss4_Curve(Pars,x)
a1=Pars(1,1);
b1=Pars(2,1);
c1=Pars(3,1);
a2=Pars(1,2);
b2=Pars(2,2);
c2=Pars(3,2);
a3=Pars(1,3);
b3=Pars(2,3);
c3=Pars(3,3);
a4=Pars(1,4);
b4=Pars(2,4);
c4=Pars(3,4);
    y = a1*exp(-((x-b1)/c1).^2) + a2*exp(-((x-b2)/c2).^2) + a3*exp(-((x-b3)/c3).^2) + a4*exp(-((x-b4)/c4).^2);
end