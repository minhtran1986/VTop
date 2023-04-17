function [E] = MatIntFnc(y,type,param)
eps = 1e-4;  %Ersatz stiffness
switch(type)
    case('SIMP')
        penal = param;
        E = eps+(1-eps)*y.^penal;
    case('SIMP-H')
        penal = param(1);
        beta = param(2);
        h = 1-exp(-beta*y)+y*exp(-beta);
        E = eps+(1-eps)*h.^penal;
    case('RAMP')
        q = param;
        E = eps+(1-eps)*y./(1+q*(1-y));
    case('RAMP-H')
        q = param(1);
        beta = param(2);
        h = 1-exp(-beta*y)+y*exp(-beta);
        E = eps+(1-eps)*h./(1+q*(1-h));
    case('LOGISTIC')
        b = param(1);
        h = param(2);
        E = exp(b*y)./(h + exp(b*y))*(h + exp(b))/exp(b);
end