function [x_fit, y_fit, fit_model, gof] = fit_gaussian(xdata, ydata, fit_points, bckgrnd_val)

syms x a b c d; % a=scale parameter, b=center of the peak, c=width of the distribution

fun = d + a*exp(-1*(x-b)^2/(c^2));  %d + (a*c/(2*pi))/((x-b)^2+(0.5*c)^2);
fun = matlabFunction(fun);

fo = fitoptions('Method','NonlinearLeastSquares');
fun_fit = fittype(fun,'independent','x','options',fo);
[M,I] = max(ydata);
x0 = [M, xdata(I), 0.01, bckgrnd_val]; % [scale param, peak location, c]
lower = [0.02*M, xdata(I)-0.1, 0.005, bckgrnd_val - 500];
upper = [0.1*M, xdata(I)+0.1, 0.05, bckgrnd_val + 500];

[fit_model, gof] = fit(xdata,ydata,fun_fit,'start',x0,'Lower',lower,'Upper',upper);
x_fit = linspace(xdata(1),xdata(end),fit_points);
y_fit = feval(fit_model,x_fit);
figure
plot(x_fit,y_fit)
hold on 
plot(xdata,ydata)