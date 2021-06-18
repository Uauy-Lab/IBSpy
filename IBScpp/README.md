This is the C / C ++ version of IBSpy, the original repository: https://github.com/Uauy-Lab/IBSpy

## How to install

requirement : Installing IBSCPP requires your GCC version to support C ++ 17 features.

The installation process is very simple, go to the project directory and execute the following command.

```shell
mkdir -p build 
cd build 
cmake -DCMAKE_BUILD_TYPE=Release ../ 
make
```

Then you can find IBScpp in the `build` folder.

## helper
```
$ ./IBScpp -h
cpp version of IBSpy,original repository:https://github.com/Uauy-Lab/IBSpy
Usage:
  IBScpp [OPTION...]

  -k, --kmers_db arg     kmers_with_strand path(required)
  -r, --reference arg    reference file path(required)
  -p, --threads arg      threads(default=1) (default: 1)
  -s, --kmer_size arg    kmer size(default=31) (default: 31)
  -w, --window_size arg  window size(default=1000) (default: 1000)
  -h, --help             print help

```