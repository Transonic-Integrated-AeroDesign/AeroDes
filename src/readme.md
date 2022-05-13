# Transfer of information

Internally this is how AeroDes modules share information between core-classes.

The arrows indicate the transfer of information and normalization for certain variables.

```mermaid
classDiagram
  direction TB
  class ADvariables {
    // [var-name] : [var-type] : [var-units]  
    + rho : double : [kg/m&sup3]
    + vinf : double : [m/s]
    + amu : double : [kg/[ms]]
  }
  class ADcanarline{
    + dClcda0 : double
    + arceff : double
    + eceff : double
    + Bc0 : double
    + xacc : double
    + cac : double
    + ac : double
    
    solveLiftingLine() void
  }
  class ADprandtline{
    + em : double
    + arm : double
    + B : double
    + cxm : double
    + cam : double
    + am : double
    + rf : double
    + cx0 : double
    //
    // polar
    //
    + alr[i] &#40aka inc[i]&#41 : double
    + ald[i] &#40aka inc[i] &times degrad&#41 : double
    + cl_al[i] &#40aka cl[i]&#41 : double
    + cd_al[i] &#40aka cd[i]&#41 : double
    + cq_al[i] &#40aka cmac[i]&#41 : double
    
    solveLiftingLine() void
  }
  class ADcanareq{
    //
    // input canard variables
    //
    + dClcda0 : double : "canard lift slope"
    + arceff : double : "canard effective aspect ratio"
    + ec : double : "canard efficiency"
    + Bc : double : [m] "canard span"
    + xacc : double : [m] "aerodynamic center of canard"
    + cac : double : [m] "mean aerodynamic chord"
    + ac : double : [m&sup2] "planform area for 2 canards"
    //
    // input main-wing variables
    //
    + em : double : "oswald efficiency"
    + arm : double : "wing+fuse aspect ratio"
    + B : double : [m] "wing span"
    + cxm : double : [m] "maximum chord/fuselage"
    + cam : double : [m] "average aerodynamic chord"
    + am : double : [m&sup2] "wing+fuse planform area"
    + rf : double : [m] "fuselage radius"
    + lf : double : [m] "length of fuselage"
    
    linearModel() void
    nonlinearModel() void
  }
  
  ADvariables "1" --o "2" ADprandtline
  ADvariables "1" --o "2" ADcanarline
  ADcanarline ..o "3" ADcanareq : ADcanarline &#8594 ADcanareq,\n dClcda0 &#8594 dClcda0, \n arceff &#8594 arceff, \n eceff &#8594 ec, \n Bc0 &#8594 Bc, \n xac &times (Bc0/2) &#8594 xacc, \n cac &times (Bc0/2) &#8594 cac, \n 2 &times ac(Bc0&sup2/4) &#8594 ac
  ADvariables --o "3" ADcanareq
  ADprandtline --o "3" ADcanareq : ADprandtline &#8594 ADcanareq, \n em &#8594 em, \n arm &#8594 armeff, \n B &#8594 B, \n cxm &times (B/2) &#8594 cxm, \n am &times (B&sup2/4), \n cam &times (B/2) &#8594 cam, \n rf &times (B/2) &#8594 rf, \n cx0 &#8594 lf
```
