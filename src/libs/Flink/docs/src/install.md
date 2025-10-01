As mentioned before, almost all of my Fortran codes depend on the *Flink* library. So, if you want to try my codes, perhaps you have to compile and install *Flink* at first. Do not be nervous. It is quite easy.

## Download

At first, you have to download the latest source codes from the official repository of [*Flink*](https://github.com/huangli712/Flink). If every thing is OK, you will have a compressed file, such as `flink.zip` or `flink.tar.gz`.

## Uncompress

And then you need to decompress it with `unzip` or `tar`.

```shell
$ unzip flink.zip
```

or

```shell
$ tar xvfz flink.tar.gz
```

## Build Library

The next step is to compile *Flink*. Please execute the following three commands in your terminal.

```shell
$ cd flink/build
$ editor make.inc
$ make
```

Once the compilation is finished (It takes about a few seconds usually), we can find `libflink.a` and a lot of `*.mod` files under the `flink/src` directory. That is all.

!!! note

    In order to compile *Flink*, you have to ensure that a MPI-enabled Fortran compiler (which is compatibe with Fortran 2018 standard), and a numerical library which implements the BLAS and LAPACK interfaces (such as Intel MKL), are installed and configured correctly in your system. You can setup your compiling environment via editing the `build/make.inc` file.

    See [make.inc](appendix/make.md) for more details.

## Build Documentation

This step is optional. To build the documentation for the *Flink* library, please execute the following commands:

```shell
$ cd flink/docs
$ julia make.jl
```

Then after a few seconds, the documentation will be ready in the *iqist/docs/build* directory. Of course, an online documentation is available in:

```text
https://huangli712.github.io/projects/flink/index.html
```

!!! tip

    To generate the documentation, the Julia interpreter must be available in your system. See [Julia language](https://julialang.org/) for more details.
