"""
    check the version of julia environment
"""
function check_version()
    if VERSION < v"1.0-"
        println("Please use julia v1.5.0+")
        exit(-1)
    else
        println("Well, julia environment is good")
    end
end
