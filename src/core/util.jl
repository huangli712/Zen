"""
    check_version()

Check the version of julia runtime environment
"""
function check_version()
    if VERSION < v"1.5-"
        error("Please use julia v1.5.0+")
    end
end

function check_home()
    try
        home = ENV["ZEN_HOME"]
    catch e
        if isa(e, KeyError)
            error("Please setup the envirnoment variable ZEN_HOME correctly")
        end
    end
    home
end
