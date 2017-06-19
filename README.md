# The `pst-ode` PSTricks package

© 2012--`\today` Alexander Grahn

https://github.com/agrahn/pst-ode

This package defines `\pstODEsolve` for solving initial value problems for sets of Ordinary Differential Equations (ODE) using the Runge-Kutta-Fehlberg (**`RKF45`**) method with automatic step size adjustment.

`\pstODEsolve[<options>]{<result>}{<output format>}{`*t*<sub>0</sub>`}{`*t*<sub>e</sub>`}{`*N*`}{`**x**<sub>0</sub>`}{`**f**`(`*t*`, `**x**`)}`

The state vectors **x** found at *N* equally spaced output points between *t*<sub>0</sub> and *t*<sub>e</sub>  are stored in the PostScript object `<result>`,  formatted according to the second argument `<output format>`, as a long list of values. `<output format>` lists the quantities to be stored in `<result>`. The user can select from the elements of **x** and the integration parameter *t*.

The right-hand side **f**`(`*t*`, `**x**`)` of the ODE system can be input in algebraic notation, if desired. RPN (Postfix) notation of PostScript can as well be used.

`<result>` can be directly used as the `<data>` argument of `\listplot{<data>}` (package `pst-plot`) or `\listplotThreeD{<data>}` (package `pst-3dplot`).

This material is subject to the [LaTeX Project Public License](LICENSE).
