using Plots

f(x) = (x.-2).^9
g(x) = (x.^9)-(18*x.^8)+(144*x.^7)-(672*x.^6)+(2016*x.^5)-(4032*x.^4)+(5376*x.^3)-(4608*x.^2)+(2304*x).-512

xspan = Array(1.920:0.001:2.080)
plot(xspan, f.(xspan), label = "f(x) = (x-2)^9")
plot!(xspan, g.(xspan), label = "g(x)")