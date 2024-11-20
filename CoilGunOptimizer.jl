using SNOW
using OrdinaryDiffEq
using Plots

#global variables
barrel_radius = 1 #cm
wire_radius = 0.06438 / 2 #cm
unit_resistance = 16.14 #milliohms/ft
unit_resistance = unit_resistance * (1/1000) * (1/12) * (1/2.54) #ohms/cm
EMF_constant = 0.1
mass = 0.01 #kg
friction_coefficient = 0.03 
max_voltage = 12

"""
voltage(max_v, cutoff_time, curr_time, current, diode_resistance, forward_voltage)
This function mimics an electric switch, providing max_v voltage until cutoff_time.

max_v is the voltage returned if curr_time is less than the cutoff_time (0 is returned otherwise)
cutoff_time is the instant when the switch is turned off
curr_time is the time at which you want to evaluate the voltage
"""
function voltage(max_v, cutoff_time, curr_time, current, diode_resistance, forward_voltage)
    if curr_time < 0.0
        return 0.0001
    elseif curr_time >= 0 && curr_time < cutoff_time
        return max_v
    elseif curr_time >= cutoff_time && diode_resistance*current > forward_voltage
        return -forward_voltage
    else
        return -diode_resistance*current #check that the diode resistance doesn't change
    end
end

#uses
function coil(min_radius, coil_length, wire_radius, wire_length, unit_resistance)
    resistance = unit_resistance * wire_length
    turns = 0
    layers = 0
    r_add = 0.0
    wire_left = wire_length
    layer_turns = round(coil_length / (2*wire_radius))
    layer_radius = min_radius
    layer_length = layer_turns * sqrt((layer_radius^2) + (wire_radius^2))
    while wire_left >= layer_length && layer_turns > 0
        layers = layers + 1
        r_add = r_add + (layer_radius * layer_turns)
        turns = turns + layer_turns
        wire_left = wire_left - layer_length
        layer_radius = layer_radius + (wire_radius*sqrt(3))
        layer_turns = layer_turns - 1
        layer_length = layer_turns * sqrt((layer_radius^2) + (layer_length^2))
    end
    average_radius = r_add / turns
    depth = layer_radius - min_radius
    inductance = ((average_radius^2)*(turns^2)) / ((19*average_radius)+(29*coil_length)+(32*depth))
    return resistance, inductance
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
    diode_resistance = 5 #ohms
    forward_voltage = 0.51 #volts
    #calculate the voltage using the voltage function
    V_a = voltage(max_v, cutoff_time, t, u[1], diode_resistance, forward_voltage)

    #this is the juicy part... the EOM (du is the derivative and u is the state vector)
    du[1] = -(R/L_a)*u[1] + sign(u[2])*(K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = -sign(u[2])*(K_e/m)*u[1] - sign(u[3])*mu_k*g
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
    coil_length = x[2] #cm
    wire_length = x[3] #cm
    resistance, inductance = coil(barrel_radius, coil_length, wire_radius, wire_length, unit_resistance)
    inductance = inductance / (10^6) #convert uH to H
    cutoff_time = x[1]
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, max_voltage, cutoff_time] #put the variables into the vector

    #println[p]

    t = (0.0, 10.0) #time 
    u0 = [0.0, -0.05, 0.0] #initial conditions [current, position, velocity]
    
    #objective (velocity when x = pos)
    pos = 0.1
    
    #set up end condition for simulation (end simulation when the mass reaches pos)
    condition(u,t,integrator) = u[2] - pos # Is zero when u[2] = pos
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
    smallest_coil = 0.08
    largest_coil = (2*0.06438*10) / (sqrt((1.5^2)+(0.06438^2)))
    x0 = [0.0, 0.5, 100]  # starting point
    lx = [0.0, smallest_coil, 5]  # lower bounds on x
    ux = [1.0, largest_coil, 20]  # upper bounds on x
    ng = 1  # number of constraints
    lg = [0.0]  # lower bounds on g
    ug = [1.0]  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(ob!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)

    return xopt
end

function plot_a_sol(cutoff_time, coil_length, wire_length)
    #set up variables
    resistance, inductance = coil(barrel_radius, coil_length, wire_radius, wire_length, unit_resistance)
    inductance = inductance / (10^6) #convert uH to H
    p = [inductance, resistance, EMF_constant, mass, friction_coefficient, max_voltage, cutoff_time] #put the variables into the vector

    t = (0.0, 0.1) #time 
    u0 = [0.0, -0.05, 0.0] #initial conditions [current, position, velocity]

    #objective (velocity when x = pos)
    pos = 0.1

    #set up end condition for simulation (end simulation when the mass reaches pos)
    #condition(u,t,integrator) = u[2] - pos # Is zero when u[2] = pos
    #affect!(integrator) = terminate!(integrator)
    #cb = ContinuousCallback(condition, affect!)
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    #sol = solve(prob, callback = cb); #solve the problem
    sol = solve(prob)
    plot(sol, label = ["current" "position" "velocity"])
end

function plot_end_velocities()
    #set up variables
    coil_length = 1 #cm
    wire_length = 100 #cm
    resistance, inductance = coil(barrel_radius, coil_length, wire_radius, wire_length, unit_resistance)
    inductance = inductance / (10^6) #convert uH to H

    times = range(start = 0.001, stop = 0.2, step = 0.001)
    velocities = []
    for t in times
        p = [inductance, resistance, EMF_constant, mass, friction_coefficient, max_voltage, t] #put the variables into the vector

        println("p: " * string(p))

        t = (0.0, 2.0) #time 
        u0 = [0.0, -0.05, 0.0] #initial conditions [current, position, velocity]

        #objective (velocity when x = pos)
        pos = 0.1

        #set up end condition for simulation (end simulation when the mass reaches pos)
        condition(u,t,integrator) = u[2] - pos # Is zero when u[2] = pos
        affect!(integrator) = terminate!(integrator)
        cb = ContinuousCallback(condition, affect!)
        #simulate
        prob = ODEProblem(single_stage, u0, t, p) #set up the problem
        sol = solve(prob, callback = cb); #solve the problem
        # get the velocity when the position is pos
        push!(velocities, sol[length(sol)][3])
    end
    plot(times, velocities)
    plot!(display = true)
end

solution = optimize_coil_gun() #this line lets you run the file to get the output
plot_a_sol(solution[1], solution[2], solution[3])

#plot_a_sol(0.01, 1, 100)

#plot_end_velocities()