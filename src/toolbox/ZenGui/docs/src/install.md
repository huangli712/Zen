## Prerequisite

The ZenGui is a Julia package. So, please make sure that the newest version of Julia language is already installed in your system. Besides the standard library, the ZenGui also relies on the following third-party packages:

* FileIO
* Images
* CImGui
* GLFW
* ModernGL

So, please install these packages by yourself. It is an easy job by using the Julia package manager or the `Pkg` package.

## How to install

Though ZenGui is not a registried Julia package, it is quite easy to install it. Perhaps the simplest way is to download it from Github. The official repository is:
```text
https://github.com/huangli712/ZenGui
```
It should be a compresed file (`ZenGui.tar.gz`). Please uncompress it in your favourite directory, such as `/home/your_home/zengui`. Then, please setup the environment variable `ZEN_GUI`, which should point to the `/home/your_home/zengui/src` directory. In MacOS, you can accomplish this by execute the following command in the terminal:
```shell
$ export ZEN_GUI=/home/your_home/zengui/src
```
Of course, you can put this command into the `.profile` or `.bashrc` file, such that the OS can execute it automatically when it is rebooted next time.

## Documentation

The official documentation for ZenGui is at the following website:

```text
https://huangli712.github.io/projects/zengui/index.html
```

Sometimes, users might want to generate the documentation by themselves. Please go to the directory `zengui/docs`, and execute the following command in the terminal:

```shell
$ julia make.jl
```

If everything is OK, you can find the generated documentation at the folder `zengui/docs/build`. Please open the `index.html` file with your favourite web browser.
