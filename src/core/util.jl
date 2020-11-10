"""
    check_version()

Check the version of julia runtime environment
"""
function check_version()
    if VERSION < v"1.5-"
        error("Please use julia v1.5.0+")
    end
end

"""
    check_home()

Check the home directory for zen
"""
function check_home()
    ENV["ZEN_HOME"]
end

"""
    check_toml()

Check whether a case.toml file is provided
"""
function check_toml()
    nargs = length(ARGS)
    if nargs < 1
       error("Please provide a case.toml file to configure your calculation")
    else
       ARGS[1]
    end
end

"""
    welcome()

Print out the welcome messages to the screen
"""
function welcome()
    println("ZEN")
    println("")
    println("version: $version")
    println("authors: $authors")
end

"""
    goodbye()

Print the goodbye messages to the screen
"""
function goodbye() 
    println("See you later")
end
