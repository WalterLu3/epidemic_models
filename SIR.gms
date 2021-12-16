$title SIR model
option limcol = 0, limrow = 0, solprint = off;
set time /0*4000/;
set state /S,I,R/;



$onexternalInput

scalar
    objectToSolve '1 is cost -1 is continuity'/-1/
    gammaValue /0.25/
    S_start /0.7/
    I_start /0.15/
    R_start /0.15/
    S_end /0.45/
    I_end /0.05/
    R_end /0.5/
;


$offExternalInput

*initial values
scalar
*   beta /0.5/
    gamma 
    h
;

gamma = gammaValue;




h = 40/(card(time)-1);

alias(time,t);


nonnegative variable
    x(state,time)
    beta(time)
*    gamma(time)
    delta(time)
    abs_max
;

beta.up(time) = 0.6;

beta.lo(time) = 0.1;

beta.fx('0') = 0.5;
x.fx('S','0') = S_start;
x.fx('I','0') = I_start;
x.fx('R','0') = R_start;

x.fx('S','4000') = S_end;
x.fx('I','4000') = I_end;
x.fx('R','4000') = R_end;
free variable
    temp
;
equation
    S_change(t)
    I_change(t)
    R_change(t)
    abs1(t)
    abs2(t)
    absMax(t)
    obj_continuity
    obj_cost
;

S_change(t)$(ord(t) < card(t))..
    x('S',t+1) =e= x('S',t) + h* ( -beta(t) *  x('S',t) * x('I',t));
    
I_change(t)$(ord(t) < card(t))..
    x('I',t+1) =e= x('I',t) + h* ( beta(t) *  x('S',t) * x('I',t) - gamma *  x('I',t));
    
R_change(t)$(ord(t) < card(t))..
    x('R',t+1) =e= x('R',t) + h* ( gamma*  x('I',t));
   
abs1(t)$(ord(t) < card(t))..
    delta(t+1) =g= beta(t+1) - beta(t);
    
abs2(t)$(ord(t) < card(t))..
    delta(t+1) =g= beta(t) - beta(t+1) ;

absMax(t)..
    abs_max =g= delta(t);
 
obj_continuity..
*    temp =g= sum( t,x('I',t) );
    temp =g= abs_max;
*    temp =g= 0;
    
obj_cost..
*    temp =g= sum( t,x('I',t) );
    temp =e= sum(t,beta(t));
*    temp =g= 0;
    

model SIR_continuity
    /
    S_change
    I_change
    R_change
    abs1
    abs2
    absMax
    obj_continuity/;
    
model SIR_cost
    /
    S_change
    I_change
    R_change
    obj_cost/;

if ( objectToSolve=1,
    solve SIR_cost using nlp maximizing temp;
else
    solve SIR_continuity using nlp minimizing temp;
);

set header /sirDetail,b/;

$onexternaloutput

parameter
result(state,t,header);

$offexternaloutput


result(state,t,'sirDetail') = x.l(state,t);
result(state,t,'b') = beta.l(t);

parameter
    endingValue(state)    
;

endingValue(state) = x.l(state,'4000');
display beta.l;
