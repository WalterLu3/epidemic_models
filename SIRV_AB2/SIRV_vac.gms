$title SIRV AB2 with Vaccination Rate as Decision Variable
option limcol = 0, limrow = 0, solprint = off;
*option nlp = baron;
*define parameter
parameter
    beta infection rate /0.5/
    alpha recovery rate /0.25/
*    u   vaccination rate /0.1/
    lambda immunity decay rate /0.05/
;

scalar
    endTime /40/
    h time interval between two mesh points
;

set
    time mesh points /t0*t8000/
    timeLessLast(time);
;
timeLessLast(time) = yes;


h = endTime/(card(time)-1);

alias (time,t);

positive variable
    S(t) Susceptible
    I(t) Infected
    R(t) Removed
    V(t) Vaccinated
    u(t) vaccinated rate
    D(t) absolute difference
    absVal
;

free variable
    objective
    
;

* bounds

u.up(t) = 0.2;

u.lo(t) = 0;
* initial value

S.fx('t0') = 0.95;
I.fx('t0') = 0.05;
R.fx('t0') = 0;
V.fx('t0') = 0;

u.fx('t0') = 0.05;


S.fx('t8000') = 0.3;
V.fx('t8000') = 0.5;

*V.lo('t4000') = 0.4;

V.lo(t)$(ord(t) > 4000) = 0.4;

*V.lo(t)$(ord(t) > 1001) = 0.2;
    
* model formulation

*------------------------ab2 method-------------------------

equation
    AB2_S(t)
    AB2_I(t)
    AB2_R(t)
    AB2_V(t)
    Euler_S
    Euler_I
    Euler_R
    Euler_V
    vaccination_change1(t)
    vaccination_change2(t)
    maxRate1(t)
    obj
;

AB2_S(t)$(ord(t) > 1 and ord(t) < card(t))..
    S(t+1) =e= S(t) + h* ( (3/2)*(-beta*S(t)*I(t) - u(t)*S(t) + lambda*R(t) + lambda*V(t))
                    - (1/2) * (-beta*S(t-1)*I(t-1) - u(t)*S(t-1) + lambda*R(t-1) + lambda*V(t-1)));
   
AB2_I(t)$(ord(t) > 1 and ord(t) < card(t))..
    I(t+1) =e= I(t) + h* ( (3/2)*( beta*S(t)*I(t) - alpha*I(t) )
                    - (1/2) * ( beta*S(t-1)*I(t-1) - alpha*I(t-1)));
                   
AB2_R(t)$(ord(t) > 1 and ord(t) < card(t))..
    R(t+1) =e= R(t) + h* ( (3/2)*( alpha*I(t) - lambda*R(t))
                    - (1/2) * ( alpha*I(t-1) - lambda*R(t-1)));
                   
AB2_V(t)$(ord(t) > 1 and ord(t) < card(t))..
    V(t+1) =e= V(t) + h* ( (3/2)*(u(t)*S(t) - lambda*V(t))
                    - (1/2) * (u(t)*S(t-1) - lambda*V(t-1)));
                    
Euler_S..
    S('t1') =e= S('t0') + h * (-beta*S('t0')*I('t0') - u('t0')*S('t0') + lambda*R('t0') + lambda*V('t0'));

Euler_I..
    I('t1') =e= I('t0') + h * ( beta*S('t0')*I('t0') - alpha*I('t0') );

Euler_R..
    R('t1') =e= R('t0') + h * ( alpha*I('t0') - lambda*R('t0') );

Euler_V..
    V('t1') =e= V('t0') + h * ( u('t0')*S('t0') - lambda*V('t0') );



vaccination_change1(t)$(ord(t) > 1)..
    D(t) =g= u(t)*S(t) -u(t-1)*S(t-1);


vaccination_change2(t)$(ord(t) > 1)..
    D(t) =g= -u(t)*S(t) + u(t-1)*S(t-1);
    
$ontext
vaccination_change1(t)$(ord(t) > 1)..
    D(t) =g= u(t) -u(t-1);


vaccination_change2(t)$(ord(t) > 1)..
    D(t) =g= -u(t) + u(t-1);
$offtext
maxRate1(t)..
    absVal =g= D(t);

obj..
    objective =e= absVal;

*obj..
*    objective =e= sum(t,u(t)*s(t));
    
model SIRV_AB2 /
    AB2_S
    AB2_I
    AB2_R
    AB2_V
    Euler_S
    Euler_I
    Euler_R
    Euler_V
    vaccination_change1
    vaccination_change2
    maxRate1
    obj/;


solve SIRV_AB2 using nlp minimizing objective;

parameter
    endS
    endI
    endR
    endV
;

endS = S.l('t8000');
endI = I.l('t8000');
endR = R.l('t8000');
endV = V.l('t8000');


display u.l;

execute_unload "vaccination_rate_02.gdx", S,I,R,V,u;