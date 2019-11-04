%%THE PROGRAM GIVES A SOLUTION FOR ONE DIMENSIONAL HEAT TRANSFER THROUGH
%%ANY CASE WITH A CONSTANT HEAT FLUX BOUNDARY CONDITION ON BOTH THE
%%BOUNDARIES IF THE OBJECT IS SYMMETRICAL AND THE CONDITIONS ARE
%%SYMMETRICAL USING THE EXPLICIT SCHEME OF TRANSIENT FINITE VOLUME METHOD.
clear all
clc
%-------------------Input Variables----------------------%
L= input('enter the length of the domain in m - ');
rho= input('density of the material kg/m^3 - ');
c= input ('specific heat capacity of the material in J/kg-K - ');
k= input('thermal conductivity of the material W/m-K - ');
alpha= 1/(rho*c);

%---------------------grid generation-----------------------%
N= input('enter the number of node points - ');
deltax = (L/N);                                             % the distance between each node point
x= (deltax/2):deltax:L;                                     % the values of the position of each node of the cells are found out for bettr calculation of the distance between them
x1= 0:deltax:L;

%----------------critical time determination----------------%
delta_t_critical = ((1/alpha)*(deltax^2))/(2*k);

%--------------Time step input----------------%
t = input('enter the time for which the temperature profile should be calculated - ');
fprintf('the critical time step value is %f\n',delta_t_critical);                  % for explicit condition to converge, there should be no negative coefficients
deltat = input('enter a time step value less than critical time step value - ');
M=(t/deltat)+1;                                                                      % the number of time steps needed for the calculation is determined and 1 is added with the intrest of having a continious array for temperature.
for a=1:M
    time(a)= (a-1)*deltat;
end

%----------------Initial conditions--------------------%
tini1 =input('enter the initial temperature at x=L/2 - ');
tini2 =input('enter the initial temperature at x=L - ');
for i=1:N
    j=1;
    if x(i)<= (L/2)
        %for these nodes the temperature remains constant till x=L/2
       temp(i,j)= tini1;
    else
        % the gradient initial condition from L/2 is incorporated by
        % interpolation.
       temp(i,j)= tini1+(((x(i)-L/2)*(tini2-tini1))/(L-(L/2)));
    end
end

%---------------boundary1 condition------------------%
qin= input('enter the heat input at L in W/m^2 - ');


%-------------solving---------------%
for j=2:M
    for i= 1:N
        %the boundary at L=0, the node point 1 has the same
        %characteristics of the domain on the left side due to symmetry.
        if x(i)== (deltax/2)
            temp(i,j)= alpha*(deltat/deltax)*((k*temp(i+1,j-1)/(x(i+1)-x(i)))+(k*temp(i,j-1)/(x(i+1)-x(i)))+(((inv(alpha)*deltax/deltat)-(k/(x(i+1)-x(i)))-(k/(x(i+1)-x(i))))*temp(i,j-1)));
        %heat flux boundary condition qin=qout at the boundary x=L
        elseif x(i) == (L-(deltax/2))
            temp(i,j)= alpha*(deltat/deltax)*(((k*temp(i,j-1)+(qin*(deltax/2)/k))/(L-x(i)))+(k*temp(i-1,j-1)/(x(i)-x(i-1)))+(((inv(alpha)*deltax/deltat)-(k/(x(i)-x(i-1)))-(k/(L-x(i))))*temp(i,j-1)));
        else
            %solving for interior nodes
            temp(i,j)= alpha*(deltat/deltax)*((k*temp(i+1,j-1)/(x(i+1)-x(i)))+(k*temp(i-1,j-1)/(x(i)-x(i-1)))+(((inv(alpha)*deltax/deltat)-(k/(x(i+1)-x(i)))-(k/(x(i)-x(i-1))))*temp(i,j-1)));
        end
    end
end

%------------------Post Processing-------------------%
plot(time,temp);
title('Temperature vs Time');
grid on;
xlabel('Time in seconds');
ylabel('Temperature in Kelvin');
figure;
plot(x,temp(:,end))
grid on;
title('Temperature vs Length');
xlabel('Length in meters');
ylabel('Temperature in Kelvin');
