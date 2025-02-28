function [x_fit, y_fit, fit_model, gof, output] = fit_sine(xdata, ydata, fit_points)

syms x a b c d; % a=scale parameter, b=freq, c=offset angle, d=DC value

fun = a*sin(b*x+c)+d;
fun = matlabFunction(fun);

fo = fitoptions('Method','NonlinearLeastSquares','TolFun',1e-6,'TolX',1e-6); %-12
fun_fit = fittype(fun,'independent','x','options',fo);
[Mx,~] = max(ydata);
[Mn,~] = min(ydata);

[fit_model_temp] = fit(xdata, ydata-mean(ydata), 'sin1');
x0 = [fit_model_temp.a1, fit_model_temp.b1, fit_model_temp.c1, mean(ydata)]; % [scale param, period, offset_angle, dc value]

% x0 = [(Mx-Mn)/2, pi, 0, mean(ydata)]; % [scale param, period, offset_angle, dc value]
lower = [0.1*(Mx-Mn)/2, 0.75*pi,   -pi, min(ydata)];
upper = [5*(Mx-Mn)/2, 1.5*pi, pi, Mx];

[fit_model, gof, output] = fit(xdata,ydata,fun_fit,'start',x0,'Lower',lower,'Upper',upper);
x_fit = linspace(xdata(1),xdata(end),fit_points);
y_fit = feval(fit_model,x_fit);