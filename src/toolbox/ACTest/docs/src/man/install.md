## Prerequisites

We chose the Julia language to develop the ACTest toolkit. Julia is an interpreted language. Hence, to run the ACTest toolkit, the latest version of the Julia interpreter must be installed on the target system. Actually, Julia's version number should not be less than 1.6. The core features of the ACTest toolkit rely only on Julia's standard library. However, if one needs to invoke the ACFlow and MiniPole toolkits for analytic continuation calculations, these toolkits must be available in the system. As for how to install and configure the ACFlow and MiniPole toolkit, please refer to the relevant papers. Additionally, if one wants to use the built-in script of ACTest to visualize the calculated results (i.e., the spectral functions), support from the CairoMakie.jl package is necessary. CairoMakie.jl is a 2D plotting program developed in Julia. It can be installed via Julia's package manager:

```shell
julia> ]
(@v1.11) pkg> add CairoMakie
```

## Main Program

The official repository of the ACTest toolkit is hosted on github. Its URL is as follows:

```text
https://github.com/huangli712/ACTest
```

Since the ACTest toolkit has not yet been registered as a regular Julia package, please type the following commands in Julia's REPL to install it:

```julia
julia> using Pkg
julia> Pkg.add("https://github.com/huangli712/ACTest")
```

If this method fails due to an unstable network connection, an offline method should be adopted. First of all, please download the compressed package for the ACTest toolkit from github. Usually it is called actest.tar.gz or actest.zip. Second, please copy it to your favorite directory (such as /home/your\_home), and then enter the following command in the terminal to decompress it:

```shell
$ tar xvfz actest.tar.gz
```

or

```shell
$ unzip actest.zip
```

Finally, we just assume that the ACTest toolkit is placed in the directory:

```text
/home/your_home/actest
```

Now we have to manually setup the environment variable **ACTEST_HOME**. It should be associated with the above directory:

```shell
export ACTEST_HOME=/home/your_home/actest/src
```

These steps just ensure that the Julia interpreter can find and import the ACTest toolkit correctly. Now it should work as expected.

## Documentation

The ACTest toolkit ships with detailed documentation, including the user's manual and application programming interface. They are developed with the Markdown language and the Documenter.jl package. So, please make sure that the latest version of the Documenter.jl package is correctly installed (please visit the website https://github.com/JuliaDocs/Documenter.jl for more details). If everything is OK, users can generate the documentation by themselves. Please enter the following commands in the terminal:

```shell
$ pwd
/home/your_home/actest
$ cd docs
$ julia make.jl
```

After a few seconds, the documentation is built and saved in the actest/docs/build directory. The entry of the documentation is actest/docs/build/index.html. You can open it with any web browser.

!!! tip

    If the ACTest toolkit was intalled in the offline mode, please make sure that the environment variable **ACTEST_HOME** is correctly setup, or else you can not launch the scripts in actest/util folder or build the documentation. A good place to setup the environment variable **ACTEST_HOME** is in the **.bashrc** (for Linux system) or **.profile** (for MacOS system) file.
