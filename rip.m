figure(1);
% f1 = -(3.8462*x1+0.2166*x3^2*sin(x2)-19.5258*cos(x2)*sin(x2)+0.00011352*x3*cos(x2)+0.1083*x1*x3*sin(x2)-0.1598*x1^2*cos(x2)^2*sin(x2)+0.1083*x1*x3*cos(x2)*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
% f2 = x3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
% f3 = -(0.000524*x3-90.1468*sin(x2)-9.7629*sin(x2)^2+5.6766*x1*cos(x2)+0.00005676*x3*sin(x2)-0.0799*x1^2*cos(x2)*sin(x2)^2-0.7379*x1^2*cos(x2)*sin(x2)+0.3197*x3^2*cos(x2)*sin(x2)+0.0160*x1*x3*cos(x2)*sin(x2)+0.1598*x1*x3*cos(x2)^2*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
% 
% g1 = 76.9231/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
% g2 = 0;
% g3 = 113.5308*cos(x2)/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
% 
% f = [f1;f2;f3];
% g = [g1;g2;g3];

u=0;
x1=0;x2=pi/4;x3=0;

for i = 0:10000
    f1 = -(3.8462*x1+0.2166*x3^2*sin(x2)-19.5258*cos(x2)*sin(x2)+0.00011352*x3*cos(x2)+0.1083*x1*x3*sin(x2)-0.1598*x1^2*cos(x2)^2*sin(x2)+0.1083*x1*x3*cos(x2)*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
    f2 = x3;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
    f3 = -(0.000524*x3-90.1468*sin(x2)-9.7629*sin(x2)^2+5.6766*x1*cos(x2)+0.00005676*x3*sin(x2)-0.0799*x1^2*cos(x2)*sin(x2)^2-0.7379*x1^2*cos(x2)*sin(x2)+0.3197*x3^2*cos(x2)*sin(x2)+0.0160*x1*x3*cos(x2)*sin(x2)+0.1598*x1*x3*cos(x2)^2*sin(x2))/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);

    g1 = 76.9231/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);
    g2 = 0;
    g3 = 113.5308*cos(x2)/(0.1083*sin(x2)-0.3197*cos(x2)^2+1);

    f = [f1;f2;f3];
    g = [g1;g2;g3];
    
    dX = f+g*u;
    x1 = x1+0.001*dX(1);
    x2 = mod(x2+0.001*dX(2),2*pi);
    if(x2>pi) x2=-2*pi+x2;  end
    x3 = x3+0.001*dX(3);
    
    plot(i,x1,'r');
    hold on;
    if(x2>0) X2=x2;
    else X2=x2+2*pi;   end
    plot(i,X2,'b');
    plot(i,x3,'g');
end
title('Free motion');
legend('x1:Dalpha','x2:beta','x3:Dbeta');
xlabel('Time(msec)'); ylabel('Value(radian)');
