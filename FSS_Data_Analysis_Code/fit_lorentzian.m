function [x_fit, y_fit, fit_model, gof] = fit_lorentzian(xdata, ydata, fit_points)

syms x a b c; % a=scale parameter, b=center of the peak, c=width of the distribution

fun = (a*c/(2*pi))/((x-b)^2+(0.5*c)^2);
fun = matlabFunction(fun);

fo = fitoptions('Method','NonlinearLeastSquares');
fun_fit = fittype(fun,'independent','x','options',fo);
[M,I] = max(ydata);
x0 = [M, xdata(I), 0.01]; % [scale param, peak location, c]
lower = [0.02*M, xdata(I)-0.1, 0.005];
upper = [0.1*M, xdata(I)+0.1, 0.05];

[fit_model, gof] = fit(xdata,ydata,fun_fit,'start',x0,'Lower',lower,'Upper',upper);
x_fit = linspace(xdata(1),xdata(end),fit_points);
y_fit = feval(fit_model,x_fit);