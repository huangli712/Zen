It is quite easy to install `ACGui`.

## Prerequisite

Besides the Julia runtime environment, `ACGui` only needs the following packages:

* TOML
* Base64
* Dash
* ACFlow

Here, `TOML` and `Base64` are standard libraries. `Dash` is a web framework, users have to install it manually.

```julia-repl
julia> using Pkg
julia> Pkg.add("Dash")
```

Note that [`ACFlow`](https://github.com/huangli712/ACFlow) is not a registered package, so please follow its user guide to install it.

## How To Install

The users can use the `Pkg` package to install the `ACGui` app from its github repository directly:

```julia-repl
julia> using Pkg
julia> Pkg.add(url = "https://github.com/huangli712/ACGui")
```

If the installed `ACGui` app is outdated, the users can use the following commands to upgrade `ACGui`:

```julia-repl
julia> using Pkg
julia> Pkg.update("ACGui")
```

!!! info

    If the users do not want to install `ACGui`, just have a try. The simplest way is as follows:

    * Download the source codes of `ACGui` from its github repository:

    ```text
    https://github.com/huangli712/ACGui
    ```

    Then uncompress it into `/home/your_home/acgui`.

    * Plug the following code in front of the `acgui/util/acg.jl` script:

    ```julia
    push!(LOAD_PATH, "/home/your_home/acgui/src")
    ```

    or just setup the environment variable `ACGUI_HOME`:

    ```shell
    export ACGUI_HOME=/home/your_home/acgui/src
    ```

## Documentation

Finally, in order to generate the documentation, please type the following commands in the terminal:

```shell
$ pwd
/home/your_home/acgui
$ cd docs
$ julia make.jl
```

After a few seconds, the documentation is built and saved in the `acgui/docs/build` directory if everything is OK. The home page of the documentation is `acgui/docs/build/index.html`. The users can open it with any web browsers.

!!! info

    If the documentation has not been generated successfully, please try the following solutions:

    * Install the `ACFlow` and `ACGui` codes via the official `Pkg` package.
    * Or setup the environment variables `ACFLOW_HOME` and `ACGUI_HOME`. They should be associated with the `src` folders of the `ACFlow` and `ACGui` codes, respectively.
