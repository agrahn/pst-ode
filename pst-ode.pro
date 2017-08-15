%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PostScript prologue for pst-ode.tex.
% Version 0.11, 2017/08/15
%
% Alexander Grahn (C) 2012--today
%
% This program can be redistributed and/or modified under the terms
% of the LaTeX Project Public License Distributed from CTAN archives
% in directory macros/latex/base/lppl.txt.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/tx@odeDict 1 dict def
tx@odeDict begin
/ode@@dict 1 dict def
/ode@dict {ode@@dict begin} bind def
ode@dict
  %some constants for step size calculation
  /sfty 0.9 def /pgrow -0.2 def /pshrink -0.25 def
  %helper functions
  /addvect { % [1 2 3] [4 5 6] addvect => [5 7 9]
    ode@dict
      aload pop xlength1 -1 roll {xlength1 -1 roll add} forall
      xlength array astore
    end
  } bind def
  /subvect { % [1 2 3] [4 5 6] subvect => [-3 -3 -3]
    ode@dict
      aload pop xlength1 -1 roll {xlength1 -1 roll sub} forall
      xlength array astore
    end
  } bind def
  /mulvect { % [1 2 3] 4 mulvect => [4 8 12]
    ode@dict /mul cvx 2 array astore cvx forall xlength array astore end
  } bind def
  /edivvect { % [1 2 3] [4 5 6] edivvect => [0.25 0.4 0.5]
    ode@dict
      aload pop xlength1 -1 roll {xlength1 -1 roll div} forall
      xlength array astore
    end
  } bind def
  /eabsvect { % [-1 2 -3] eabsvect => [1 2 3]
    ode@dict {abs} forall xlength array astore end
  } bind def
  %/revstack { % (a) (b) (c) (d) 3 revstack => (a) (d) (c) (b)
  %  -1 2 {dup 1 sub neg roll} for
  %} bind def
  /min { 2 copy gt { exch } if pop } bind def
  /max { 2 copy lt { exch } if pop } bind def
  %coefficient table (Butcher table) of RKF45
  /a41 25 216 div def /a43 1408 2565 div def
  /a44 2197 4104 div def /a45 1 5 div neg def
  /a51 16 135 div def /a53 6656 12825 div def
  /a54 28561 56430 div def /a55 9 50 div neg def
  /a56 2 55 div def
  /b21 1 4 div def /b31 3 32 div def /b32 9 32 div def
  /b41 1932 2197 div def /b42 7200 2197 div neg def /b43 7296 2197 div def
  /b51 439 216 div def /b52 8 neg def /b53 3680 513 div def
  /b54 845 4104 div neg def /b61 8 27 div neg def /b62 2 def
  /b63 3544 2565 div neg def /b64 1859 4104 div def /b65 11 40 div neg def
end
%Runge-Kutta-Fehlberg (RKF45) method
%performs one integration step over tentative step size ddt
%[state vector x(t)] RKF45 => [x(t)] [x(t+ddt) by RKF4] errmax
/RKF45 {
  dup
  ode@dict tcur end /t exch def
  ODESET
  ode@dict
    ddt mulvect /k1 exch def
    dup k1 b21 mulvect addvect
    tcur ddt 4 div add
  end
  /t exch def
  ODESET
  ode@dict
    ddt mulvect /k2 exch def
    dup k1 b31 mulvect addvect k2 b32 mulvect addvect
    tcur ddt 3 mul 8 div add
  end
  /t exch def
  ODESET
  ode@dict
    ddt mulvect /k3 exch def
    dup k1 b41 mulvect addvect k2 b42 mulvect addvect k3 b43 mulvect addvect
    tcur ddt 12 mul 13 div add
  end
  /t exch def
  ODESET
  ode@dict
    ddt mulvect /k4 exch def
    dup k1 b51 mulvect addvect k2 b52 mulvect addvect k3 b53 mulvect addvect
    k4 b54 mulvect addvect
    tcur ddt add
  end
  /t exch def
  ODESET
  ode@dict
    ddt mulvect /k5 exch def
    dup k1 b61 mulvect addvect k2 b62 mulvect addvect k3 b63 mulvect addvect
    k4 b64 mulvect addvect k5 b65 mulvect addvect
    tcur ddt 2 div add
  end
  /t exch def
  ODESET
  ode@dict
    ddt mulvect /k6 exch def % => [x(t)]
    %fourth order solution (increment dx)
    dup dup k1 a41 mulvect k3 a43 mulvect addvect k4 a44 mulvect addvect
      k5 a45 mulvect addvect dup
      % => [x(t)] [x(t)] [x(t)] [dx by RKF4] [dx by RKF4]
    %fifth order solution (abs. error)
    k1 a51 mulvect k3 a53 mulvect addvect k4 a54 mulvect addvect
      k5 a55 mulvect addvect k6 a56 mulvect addvect subvect
      % => [x(t)] [x(t)] [x(t)] [dx by RKF4] [err]
    5 1 roll addvect 4 -2 roll % => [x(t)] [x(t+ddt) by RKF4] [err] [x(t)]
    %scaling vector for step size adjustment (Numerical Recipies)
    eabsvect k1 eabsvect addvect {1e-30 add} forall xlength array astore
    % => [x(t)] [x(t+ddt) by RKF4] [err] [xscale]
    edivvect eabsvect 1e-30 exch {max} forall %maximum rel. error
    % => [x(t)] [x(t+ddt) by RKF4] errmax
  end
} bind def
/ODEINT { % performs integration over output step [t,t+dt]
          % [ state vector x(t) ] ODEINT => [ state vector x(t+dt) ]
  %decrease overshooting step size
  ode@dict
    tout tcur sub dup abs ddt abs lt {/ddt exch def}{pop} ifelse
  end
  RKF45
  ode@tol div dup 1 gt {
    %failed step -> reduce step size
    ode@dict
      exch pop pshrink exp 0.1 max sfty mul ddt mul /ddt exch def
      ode@dict tcur ddt add tcur end eq {
        % error: step size underflow in ODEINT
        (! t=) odeprint tcur odeprint
        % print & remove previous state vector
        (, x=[) odeprint /ode@spc () def
        {ode@spc odeprint /ode@spc ( ) def odeprint} forall (]) odeprint
        true
      }{
        (-) odeprint
        false
      } ifelse
    end
    % on step size underflow, leave loop over output steps (pst-ode.tex)
    {exit} if
    ODEINT %repeat step with new ddt
  }{
    %success
    3 -1 roll pop % remove previous state vector
    ode@dict
      /tcur tcur ddt add def
      dup 0 ne {pgrow exp 5 min sfty mul ddt mul /ddt exch def}{pop} ifelse
      tcur tout sub %output step completed?
    end
    0 ge {
      (o) odeprint
      ode@dict
        /tcur tout def
        /outStepCnt outStepCnt 1 add def
        /tout tStart dt outStepCnt mul add def
      end
    }{
      (+) odeprint ODEINT %continue integration
    } ifelse
  }ifelse
} bind def
/writeresult { %write output vector to file
  /loopproc load forall
  statefile (\n) writestring
} bind def
/loopproc {
  0 cvs statefile exch writestring
  statefile ( ) writestring
} bind def
/loopproc load 0 256 string put
/assembleresult { %assembles state vector for building table of results
  {
    dup (t) eq {
      pop ode@dict tcur end
    }{
      ode@laststate exch get
    } ifelse
  } forall
} bind def
end
