usage
=====

## label propagation

### Description

An implementation of Label propagation ([Zhu and Ghahramani, 2002]).

### Installation

    xiaojin$ wget http://github.com/smly/label-propagation/tarball/0.1.1 -O lprop-0.1.1.tar.gz
    xiaojin$ tar zxvf lprop-0.1.1.tar.gz


    xiaojin$ ./autogen.sh
    xiaojin$ ./configure
    xiaojin$ make

### Input format

Take it easy.

### Sample input

    xiaojin$ cat dat/sample_data7/sample.in
    2:1 3:1
    1:1 3:1
    1:1 2:1 4:1
    3:1 5:1 8:1
    4:1 6:1
    5:1 7:1
    6:1
    4:1 9:1
    8:1


    xiaojin$ cat dat/sample_data7/sample.label
    ?
    1
    ?
    ?
    ?
    2
    ?
    ?
    2

### Sample output

    xiaojin$ ./lprop -e 1e-10 -p 9 -m 1000 -i dat/smaple_data7/sample.in \
                 -l dat/sample_data7/sample.label \
                 -r dat/sample_data7/return
    Number of nodes:              9
    Number of labeled nodes:      3
    Number of unlabeled nodes:    6
    eps:                          1e-10
    max iteration:                1000
    ........................
    done
    iteration: 47 times, error: 7.76231e-11


    xiaojin$ cat dat/sample_data7/return.result
    1
    1
    1
    2
    2
    2
    2
    2
    2


    xiaojin$ cat dat/sample_data7/return.weight
    U: 0.912762520 0.087237480
    L: 1.000000000 0.000000000
    U: 0.781906300 0.218093700
    U: 0.258481422 0.741518578
    U: 0.103392569 0.896607431
    L: 0.000000000 1.000000000
    U: 0.000000000 1.000000000
    U: 0.064620355 0.935379645
    L: 0.000000000 1.000000000

### Author

Kohei Ozaki (eowner@gmail.com)
