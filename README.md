# additive-fft
Additive Fast Fourier Transforms over Finite Fields

## To-do:
1. Implement Gao's general algorithm using Cantor's basis.
1. Implement Gao's special case algorithm
1. Implement the IFFT algorithm for Cantor's algorithm in Sage
1. Implement the IFFT algorithm for Gao's algorithm in Sage 
1. Design Gao's and Cantor's algorithms for the affine space
1. Implement Gao's and Cantor's algorithms for the affine space. I should change the current implement to general affine space implementations.

### Done:
 1. Make `S_shifts_table_generator` function faster in the pre computations in the Cantor's algorithm using the fact that $S_{i+k}(y_k) = y_i$.
 1. Write comments on the code
 1. Gao's FFT excluding the pre-computations


## Timings

### 1. $m = 14, GF(2^{256})$:
#### Cantor's FFT (average over 100 iterations)
- excludes pre-computation: 0.5657543277740479 s
- includes pre-computation: 0.7454671454429627 s
#### Gao's FFT (average over 10 iterations)
- includes pre-computation: 5.390587258338928 s

