
using HG1G2
using HG1G2conversions
using Winston

rad(x) = pi*x/180

function main()
    data = readcsv("44_Nysa.csv")
    basis = default_basis()
    X = linspace(0.0, 60.0, 100)
    T = rad(X)
    tic()
    A = fit_HG1G2(data, basis, true)
    G = a1a2a3_to_HG1G2(A)
    toc()
    println(A)
    println(G)
    Y = fitted_curve(T, G, basis)
    plot(data[:,1],data[:,2],"x",X,Y)
end
