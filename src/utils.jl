function now_string(; hour=false, minute=false)
    fmt = "YYYY-mm-dd"
    if minute
        fmt *= " HH:MM"
    elseif hour
        fmt *= " HH"
    end
    return Dates.format(now(), fmt)
end
