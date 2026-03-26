#############导入常用同元库函数#####
using TyBase #基础库
using TyMath #数学库
using TyPlot #图形库
using TyStatistics #统计库
using TyCurveFitting #曲线拟合库
using TySymbolicMath #符号计算库
using TyOptimization #优化库
using TyGlobalOptimization #全局优化库
using TyMachineLearning #机器学习库
####################################

# Use a fixed-step RK4 integrator so the script does not rely on extra ODE packages.
function rk4_system(fun, tspan, x0, n::Int=200)
    t0, t1 = tspan
    h = (t1 - t0) / n
    t = collect(range(t0, t1, length=n + 1))
    y = zeros(n + 1, length(x0))
    y[1, :] = x0

    for i in 1:n
        ti = t[i]
        xi = vec(y[i, :])
        k1 = fun(ti, xi)
        k2 = fun(ti + h / 2, xi + (h / 2) * k1)
        k3 = fun(ti + h / 2, xi + (h / 2) * k2)
        k4 = fun(ti + h, xi + h * k3)
        y[i + 1, :] = xi + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
    end

    return t, y
end

function shooting(p, q, f, tspan, x0f; n::Int=200)
    ga, gb = x0f
    f1 = (t, x) -> [x[2], -q(t) * x[1] - p(t) * x[2]]
    f2 = (t, x) -> begin
        base = f1(t, x)
        [base[1], base[2] + f(t)]
    end

    t, y1 = rk4_system(f1, tspan, [1.0, 0.0], n)
    _, y2 = rk4_system(f1, tspan, [0.0, 1.0], n)
    _, yp = rk4_system(f2, tspan, [0.0, 0.0], n)

    m = (gb - ga * y1[end, 1] - yp[end, 1]) / y2[end, 1]
    t, y = rk4_system(f2, tspan, [ga, m], n)
    return t, y
end

p(t) = 1.0
q(t) = 1.0
f(t) = -1.0
x0 = [0.0, 0.0]
t, y = shooting(p, q, f, (0.0, 1.0), x0; n=200)

x1 = collect(0.0:1.0 / 56.0:1.0)
y1 = cos.(x1) .+ (1.0 - cos(1.0)) / sin(1.0) .* sin.(x1) .- 1.0

hold("on")
plot(x1, y1, "or")
plot(t, y[:, 1], "-b")
