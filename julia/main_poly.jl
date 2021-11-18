
using Printf

@inline function horner(x, p)
    hi = p[end]
    
    for i in length(p)-1:-1:1
        pi = p[i]
        prod = hi * x
        hi = pi + prod
    end
    return hi
end


@inline function two_sum(a, b)
    hi = a + b
    a1 = hi - b
    b1 = hi - a1
    lo = (a - a1) + (b - b1)
    return hi, lo
end

# return x, y, floats,  with x + y = a * b
@inline function two_product_fma(a, b)
    x = a * b
    y = fma(a, b, -x)
    x, y
end

@inline function exthorner(x, p)
    hi, lo = p[end], zero(x)
    
    for i in length(p)-1:-1:1
        pi = p[i]
        # prod,err1 = two_product_fma(hi, x)
        prod = hi * x
        err1 = fma(hi, x, -prod)
        
        hi,err2 = two_sum(pi, prod)
        
        lo = fma(lo, x, err1 + err2)
    end
    return hi, lo
end

# untested
@inline function compensated_horner(ps, x) 
    n, T = length(ps), eltype(ps)
    aᵢ = ps[end]
    sᵢ = aᵢ * _one(x)
    c = zero(T) * _one(x)
    for i in n-1:-1:1
	aᵢ = ps[i]
        pᵢ, πᵢ = two_product_fma(sᵢ, x)
	sᵢ, σᵢ = two_sum(pᵢ, aᵢ)
        c = fma(c, x, πᵢ + σᵢ)
    end
    sᵢ + c
end



function test(coeff, t::Float64)
    x = evalpoly(t, coeff)
    # x = horner(t, coeff)
    
    y = exthorner(t, coeff)
    
    delta = abs(x - y[1])
    nu = delta / eps(x)
    @printf("t: %.20f horner: %.20f exthorner: (%.20f, %.20f) delta: %.20f  ulp: %.20f\n",
                t, x, y[1], y[2], delta, nu) 
end




function solveQuad(dcoeff)
    #    d:      1.4899964868068868E19
    #    sqrt(d): 3.8600472624138775E9
    #    q: 4.919281402925132E9
    #    r2: 0.16051708637880857
    a = dcoeff[3]
    b = dcoeff[2]
    c = dcoeff[1]

    d = b * b - 4.0 * a * c
    println("d: ", d)
    d = sqrt(d)
    println("sqrt(d): ", d)

    if (b < 0.0)
        d = -d
    end

    q = -0.5 * (b + d)
    println("q: ", q)
    
    r = q / a   # first root (bigger in magnitude)
    println("r: ", r)
    
    if (q != 0.0)
        r2 = c / q
        println("r2: ", r2)
        r,r2
    end
    r
end


function solveQuad_HA(dcoeff)
    #    d:      1.4899964868068868E19
    #    sqrt(d): 3.8600472624138775E9
    #    q: 4.919281402925132E9
    #    r2: 0.16051708637880857
    a = dcoeff[3]
    b = dcoeff[2]
    c = dcoeff[1]

    #d = b * b - 4.0 * a * c
    b2,eb2 = two_product_fma(b, b)
    ac,eac = two_product_fma(-4.0 * a, c)

    d,ed = two_sum(b2, ac)
    ed += eb2 + eac
    println("d: ", d)
    println("ed: ", ed)
    
    d = sqrt(d)
    println("sqrt(d): ", d)

    if (b < 0.0)
        d = -d
    end

    q = -0.5 * (b + d)
    println("q: ", q)
    
    r = q / a   # first root (bigger in magnitude)
    println("r: ", r)
    
    if (q != 0.0)
        r2 = c / q
        println("r2: ", r2)
        return r,r2
    end
    return r
end


function main()

    N = 100
    
    coeff = [0, 789628717.875, -2989257771.71819305419921875, 2199629053.84319305419921875]
    dcoeff = [789628717.875, -5978515543.4363861083984375, 6598887161.52957916259765625]

    println("coeff: ", coeff)
    println("dcoeff: ", dcoeff)
    
    roots = solveQuad(dcoeff)
    println("roots:    ", roots)
    
    roots_HA = solveQuad_HA(dcoeff)
    println("roots_HA: ", roots_HA)
    
    for i in 1:length(roots_HA)
        root = roots_HA[i]
        @printf("root: %.20f\n", root)
        test(coeff, root)    
    end
    
    # test exthorner = 0.0
    root_exthorner_eq_0 = 0.16051708637880857
    println("deriv_poly(root) = 0.0 ?")
    test(dcoeff, root_exthorner_eq_0)    
    
    println("poly(root): ")
    test(coeff, root_exthorner_eq_0)
    
    if (false)
        # iterate from roots[2]
        real_root = 0.161
        tr = roots_HA[2]
        deltaEps = (real_root - tr) / 1000.0
        
        while (tr < real_root)
            test(dcoeff, tr)    
            tr += deltaEps
        end
    end
    
    # java roots:
    # tExtrema:   [0.7454713624448268713241459, 0.161]
    # tExtrema_d: [0.7454713624448269, 0.16051708637880857]
    
    r1 = 0.161
    test(coeff, r1)

    println("polynom [0..1]")

    for i in 1:N
        t = i * (1.0 / N)
        test(coeff, t)
    end
end

main()

