function dx = monoD(t,x) 
dx(1,1) = 10.000040783334043226204812526703*x(2) - 9.9998480764079431537538766860962*x(1) + 0.00035946010194382038704929982486647*x(3) - 0.0044888224382146901803025684785098; 
dx(2,1) = 28.00005851710011484101414680481*x(1) - 0.99999587396337119571398943662643*x(2) - 1.0000005134443199494853615760803*x(1)*x(3) - 0.0012482798041204556938055247883312; 
dx(3,1) = 1.0000200369024696556152775883675*x(1)*x(2) - 2.6672649288293541758321225643158*x(3) - 0.00027605477827286062364464669371955*x(1) + 0.0075597122773718439248114009387791; 
end