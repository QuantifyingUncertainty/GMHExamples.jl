sincosdata(t,paras,n) = data(:array,t,applynoise!(n,sincos(t,paras)))

function sincosmodel(time,modelvals,variance,paraminit...)
    p = parameters([:a,:b],paraminit...)
    n = noise(:gaussian,variance)
    d = sincosdata(time,modelvals,n)
    model(:target,p,d,n,sincos,time;name="Sine-Cosine")
end

function sincosmodel!(time,modelvals,variance,paraminit...)
    p = parameters([:a,:b],paraminit...)
    n = noise(:gaussian,variance)
    d = sincosdata(time,modelvals,n)
    model(:target!,p,d,n,sincos!,time;name="Sine-Cosine!")
end
