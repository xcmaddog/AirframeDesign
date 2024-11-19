using SNOW
using OrdinaryDiffEq
using Plots

"""
voltage(max_v, cutoff_time, curr_time)
This function mimics an electric switch, providing max_v voltage until cutoff_time.

max_v is the voltage returned if curr_time is less than the cutoff_time (0 is returned otherwise)
cutoff_time is the instant when the switch is turned off
curr_time is the time at which you want to evaluate the voltage
"""
function voltage(max_v, cutoff_time, curr_time)
    if curr_time >= cutoff_time
        return 0.0
    else
        return max_v
    end
end

#ODE for EOM of a single stage coil gun
"""
This function has the EOM for a simple, single-stage coil gun. It is used in the ODE solver
"""
function single_stage(du, u, p, t)
    #pull the variables out of the vector p, and give them more legible names
    L_a = p[1]
    R = p[2]
    K_e = p[3]#assuming K_e and K_f are the same
    m = p[4]
    mu_k = p[5]
    g = 9.81
    max_v = p[6]
    cutoff_time = p[7]
    #calculate the voltage using the voltage function
    V_a = voltage(max_v, cutoff_time, t)

    #this is the juicy part... the EOM (du is the derivative and u is the state vector)
    du[1] = -(R/L_a)*u[1] - (K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = (K_e/m)*u[1] - mu_k*g
end

"""
function ob!(g, x)
This function is used in the optimizer. 
The vector g is a vector of constraint values and x is a vector of variables the optimizer can vary to find the optimum

As it is currently written, it restricts the cutoff time to be within 0 and 2 seconds, and optimizes the cutoff time
    to find the greatest velocity at some position (pos)
"""
#objective function to optimize coil gun
#let's vary: cutoff_time
function ob!(g, x)
    #set up variables
    L_a = 1.0
    R = 1.5
    K_e = 0.1
    m = 0.01
    mu_k = 0.03
    max_v = 1500
    cutoff_time = x[1]
    p = [L_a, R, K_e, m, mu_k, max_v, cutoff_time] #put the variables into the vector
    t = (0.0, 2.0) #time 
    u0 = [0.0, -0.1, 0.0] #initial conditions [current, position, velocity]
    
    #objective (velocity when x = pos)
    pos = 0.1
    
    #set up end condition for simulation (end simulation when the mass reaches pos)
    condition(u,t,integrator) = u[2] - 0.1 # Is zero when u[2] = 0.1
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob, callback = cb); #solve the problem

    # get the velocity when the position is pos
    f = -sol[length(sol)][3] #this is negative because we are minimizing in the optimizer

    #constraints
    g[1] = x[1] #the voltage can't be turned off at a negative value nor one larger than 2s
    
    #return
    return f
end

"""
optimize_coil_gun()
This function runs, without any inputs, the optimizer for the single stage coil gun.

To read the results, look for "xstar = _____" This is the solution the optimizer found. (When to turn off the voltage)
    "fstar = _____" is the value of the objective. (The velocity of the projectile at the end of the barrel)
"""
function optimize_coil_gun()
    x0 = [0.0]  # starting point
    lx = [0.0]  # lower bounds on x
    ux = [2.0]  # upper bounds on x
    ng = 1  # number of constraints
    lg = [0.0]  # lower bounds on g
    ug = [2.0]  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(ob!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)

    return xopt[1]
end

function plot_a_sol(cutoff_time)
    #set up variables
    L_a = 0.00000534
    R = 1.5
    K_e = 0.1
    m = 0.01
    mu_k = 0.03
    max_v = 1500

    p = [L_a, R, K_e, m, mu_k, max_v, cutoff_time] #put the variables into the vector
    t = (0.0, 2.0) #time 
    u0 = [0.0, -0.1, 0.0] #initial conditions [current, position, velocity]

    #objective (velocity when x = pos)
    pos = 0.1

    #set up end condition for simulation (end simulation when the mass reaches pos)
    condition(u,t,integrator) = u[2] - pos # Is zero when u[2] = pos
    affect!(integrator) = terminate!(integrator)
    cb = ContinuousCallback(condition, affect!)
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob, callback = cb); #solve the problem
    plot(sol, label = ["current" "position" "velocity"])
end

solution = optimize_coil_gun() #this line lets you run the file to get the output

plot_a_sol(solution)