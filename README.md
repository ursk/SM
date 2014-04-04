SM
==

Simple demo for Score Matching estimation of an Energy-based (i.e. non-normalized) model. The demo is in Matlab. The example we consider is an over-complete ICA model of the form 

    s = g(Wx)
    p(s) = 1/Z exp(-E(s))

where the rows of the matrix W are filters that are applied to an image patch x, and we maximize the independence of the sources s by maximizing the likelihood under some super-gaussian energy function E. 

The demo script should be called like

    sm_ica(0, 'my_score_matching_model', .1)  

and will run for a long time, so to check progress it is advised to open a second Matlab window and run 

    load my_score_matching_model.mat; figure(1); plot(obj); figure(2); visual(dwM*V, 2, 10, 0)

after a few 100 iterations, the objective function should reduce to -45 or so, and the basis functions should start to become Gabor-like. It's possible to warmstart the estimation by passing a 1 as the first parameter to sm_ica.mat, which can be useful for on the fly tweaks to the step size, which is the 3rd parameter.