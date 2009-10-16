usage
=====

## label propagation

### Description

An implementation of Label propagation ([Zhu and Ghahramani, 2002]).

### Installation

    xiaojin$ git clone git://github.com/smly/label-propagation.git
    xiaojin$ cd label-propagation
    xiaojin$ ./configure --prefix=$HOME
    xiaojin$ make

### Input format

Take it easy.

### Sample input

    xiaojin$ cat dat/sample_data7/sample.in
    2:1 3:1
    1:1 3:1
    1:1 2:1 4:1
    3:1 5:1 8:1
    4:1 6:1 7:1
    5:1 7:1
    5:1 6:1
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

    xiaojin$ lprop -e 1e-10 -p 9 -m 1000 -i dat/smaple_data7/sample.in \
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
    U: 0.915089514 0.084910486
    L: 1.000000000 0.000000000
    U: 0.787723785 0.212276215
    U: 0.278260870 0.721739130
    U: 0.081841432 0.918158568
    L: 0.000000000 1.000000000
    U: 0.032736573 0.967263427
    U: 0.069565217 0.930434783
    L: 0.000000000 1.000000000

### Output Visualization

![graph7](http://github.com/smly/label-propagation/raw/master/dat/ss.png)

### License

BSE3 License

### Author

Kohei Ozaki (eowner@gmail.com)
