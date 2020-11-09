"""
    check_version()

Check the version of julia runtime environment
"""
function check_version()
    if VERSION < v"1.5-"
        error("Please use julia v1.5.0+")
    end
end

function check_zen_home()
end
