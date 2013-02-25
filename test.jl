
using HG1G2
using HG1G2conversions
using Winston

rad(x) = pi*x/180

function main()
    data = readcsv("44_Nysa.csv")
    data[:,1] = data[:,1] * pi/180
    basis = default_basis()
    X = linspace(0.0, max(data[:,1]), 100)
    tic()
    G = fit_HG1G2(basis, data)
    A = HG1G2_to_a1a2a3(G)
    toc()
    println(A)
    println(G)
    Y = fitted_curve(X, G, basis)
    plot(data[:,1],data[:,2],"x",X,Y)
end
