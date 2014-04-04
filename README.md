SM
==

Overcomplte ICA
=================

Simple demo for Score Matching estimation of an Energy-based (i.e. non-normalized) model. The demo is in Matlab. The example we consider is an over-complete ICA model of the form 

    s = g(Wx)
    p(s) = 1/Z exp(-E(s))

where the rows of the matrix W are filters that are applied to an image patch x, and we maximize the independence of the sources s by maximizing the likelihood under some super-gaussian energy function E. 

The demo script should be called like

    sm_ica(0, 'my_score_matching_model', .1)  

and will run for a long time, so to check progress it is advised to open a second Matlab window and run 

    load my_score_matching_model.mat; figure(1); plot(obj); figure(2); visual(dwM*V, 2, 10, 0)

after a few 100 iterations, the objective function should reduce to -45 or so, and the basis functions should start to become Gabor-like. Because of PCA whitening, the filters need to be projected back into pixel space to visualize. It's possible to warm-start the estimation by passing a 1 as the first parameter to sm_ica.mat, which can be useful for on the fly tweaks to the step size, which is the 3rd parameter.



Markov Random Field (FoE)
==========================

Convolutional energy based model, similar to the Field of Experts by Roth and Black. The demo estimates an MRF on "images" of size 24x24 pixels, which are tiled with filters of 8x8 pixels. The model is estimated with 

    cmrf(0, 'my_mrf_model', 4)

and will periodically plot output, that can be reproduced with

    load test_cmrf; figure(1); plot(obj); figure(2); visual(V, 2, 5, 0)

note that the MRF uses zero-phase whitening so the filters are estimated in pixel space and can be visualized without dewhitening. 
 
