###Target function, consistings of a 2-parameter sine-cosine function
sincos(t,paras) = hcat(sin(2*pi*paras[1]*t),cos(2*pi*paras[2]*t))

###Same target function that does not allocate memory but works in-place
function sincos!(result,timepoints,paras)
    a = 2*pi*paras[1]
    b = 2*pi*paras[2]
    @simd for i=1:length(timepoints)
        @inbounds result[i,1] = sin(a*timepoints[i])
        @inbounds result[i,2] = cos(b*timepoints[i])
    end
    result
end

