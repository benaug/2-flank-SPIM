2-flank SPIM from Augustine et al. (2018) implemented in Nimble.

https://projecteuclid.org/journals/annals-of-applied-statistics/volume-12/issue-1/Spatial-capturerecapture-with-partial-identity--An-application-to-camera/10.1214/17-AOAS1091.full

See test scripts for code comments, files used, etc. There are currently 3 versions

1) 2flank-testscript: used the full 3D data array, allows occasion effects and a variable number of operable cameras at the same site across occasions, as can happen when 1 of 2 cameras fails. Using the 3D data is slow. This version uses a rectangular state space.
2) 2flank-testscript-HabMask: Same as 1, but allows a discrete habitat mask. Space is still treated as continuous.
3) 2flank-testscript-y2D-HabMask: Same as 2, but uses 2D data array. Does not allow occasion effects or changing camera numbers, but is much faster if there are a lot of occasions. Each station has a fixed camera number and that number of cameras is either operable or not on each occasion.
