
A = rand(1_000_000_000)

function init_threads(A)
    return zeros(eltype(A), nthreads())
end

function update_R!(R, A)
    @threads for i in eachindex(A)
        @inbounds R[threadid()] += A[i]
    end
end

function gather_summands(R)
    return sum(R)
end

function threaded_sum(A)
    R = init_threads(A)
    update_R!(R, A)
    return gather_summands(R)
end

function my_sum(A)
    r = 0.
    for i in eachindex(A)
        r += A[i]
    end
    return r
end

bt = @btime threaded_sum($A)
bn = @btime my_sum(A)