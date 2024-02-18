# using Plots
using PlotlyJS

using RowEchelon
# using LinearAlgebra


P = [1.0	2.0	 3.0	1.0;
-3.0 	2.0 	5.0 	1.0;
3.141592653589793 	2.718281828459045 	-1.414213562373095 	 1.0]

A = [1.0;2.0;3.0]
B = [-3.0;	2.0;	5.0]
C = [3.141592653589793;	2.718281828459045;	-1.414213562373095]
# n = cross(C-B, A-B)
# d = dot(A, n)

d = -12.8
rP = rref(P)
a = -rP[1,4]*d
b = -rP[2,4]*d
c = -rP[3,4]*d

# n = [a,b,c]

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

# # Linear discriminant function
# p(x,y,n,m) = sin.(n*pi*x).*sin.(m*pi*y)
xs = collect((-4/0.01):(4/0.01))*0.01
ys = collect((-4/0.01):(4/0.01))*0.01
X, Y = meshgrid(xs,ys)

# a*x + b*y + c*z = -d
Z = -d .-(a/c)*X .- (b/c)*Y

layout = Layout(
    scene=attr(
        xaxis=attr(
            range=[-4,4]
        ),
        yaxis=attr(
            range=[-4,4]
        ),
        zaxis=attr(
            range=[-20,20]
        ),
    ),
    legend=attr(
        yanchor="top",
        xanchor="right",
        orientation="v"
    )
)

trace1 = surface(
        x=X,
        y=Y,
        z=Z,
        showscale=false,
    )

trace2 = scatter(
    x=P[:,1],
    y=P[:,2],
    z=P[:,3],
    mode="markers",
    marker=attr(
        size=12,
        colorscale="Viridis",   # choose a colorscale
        opacity=0.8
    ),
    type="scatter3d"
)

p = plot([trace1, trace2], layout)
savefig(p, "../figures/plane_plot.html")



