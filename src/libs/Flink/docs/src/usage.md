It is quite easy to use the *Flink* library.

## Step 1

At first, please look at the user's guide to learn how to call the *Flink*'s subroutines and functions. Notice that some subroutines are encapsulated in various modules, while the others are not.

## Step 2

Second, modify your source codes carefully.

## Step 3

Third, if your codes are using some modules of the *Flink* library, please copy the corresponding `*.mod` files to your working directory. For example, if your codes need the module `spring` which is defined in the `src/m_spring.f90` file, please copy the `spring.mod` file to your directory or add it to the INCLUDE_PATH.

## Step 4

Fourth, modify your Makefile. Notify the linker to link `libflink.a` to your program. You can use the following syntax:

```shell
$ ifort your_code.f90 -Lpath_to_libflink -lflink -o your_code.x
```
