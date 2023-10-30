# fastprimes
Fast small primes generator.

# What it is for
Generate a compressed sieve of Eratosthenes using all cores available.

* Ask for the next_prime, starting with 2.
* Test quickly if a given integer (below a certain value N) is a prime.

# Low memory usage
Compressed here means that the sieve only contains bits for integers
that are not divisible by 2, 3, 5, 7, 11 or 13. This means that of
every 30030 integers only 5760 bits are used. For example, to generate
a sieve for all integers up till 1,000,000,000,000 costs 24 GB of RAM
(23,976,023,976 bytes).

# Usage
For an example see [sum_first_n_primes.cxx](https://github.com/CarloWood/fastprimes-testsuite/blob/master/sum_first_n_primes.cxx).

Note that five other git submodules of the same author are required:
* [cwm4](https://github.com/CarloWood/cwm4)
* [cwds](https://github.com/CarloWood/cwds)
* [utils](https://github.com/CarloWood/ai-utils)
* [threadsafe](https://github.com/CarloWood/threadsafe)
* [threadpool](https://github.com/CarloWood/threadpool)

These are required for the debug code and the thread pool.

Also note that `sum_first_n_primes.cxx` is the *only* good example in the testsuite;
the rest is research that led up to this code and mostly doesn't even work.
