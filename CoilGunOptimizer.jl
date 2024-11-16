#binary search for index of lowest number exceeding some position
function find_from_sol(sol, target) 
    start = 1 
    fin = length(sol) 
    while(start <= fin)
        mid = (start + fin) ÷ 2 
        if (sol[mid][1][2] <= target)
            start = mid + 1 
        else
            fin = mid - 1
        end
    end
    return start;
end

#binary search for index of lowest number exceeding some position
function find_from_loc(loc, target)
    start = 1 
    fin = length(loc) 
    while(start <= fin)
        mid = (start + fin) ÷ 2 
        if (loc[mid] <= target)
            start = mid + 1 
        else
            fin = mid - 1
        end
    end
    return start;
end

#returns the voltage as if it was on until it was turned off at a variable time
function voltage(max_v, cutoff_time, curr_time)
    if curr_time >= cutoff_time
        return 0.0
    else
        return max_v
    end
end

#ODE for EOM of a single stage coil gun
function single_stage(du, u, p, t)
    L_a = p[1]
    R = p[2]
    K_e = p[3]#assuming K_e and K_f are the same
    m = p[4]
    mu_k = p[5]
    g = 9.81
    max_v = p[6]
    cutoff_time = p[7]
    V_a = voltage(max_v, cutoff_time, t)

    du[1] = -(R/L_a)*u[1] - (K_e/L_a)*u[3] + (V_a/L_a)
    du[2] = u[3]
    du[3] = (K_e/m)*u[1] + mu_k*g
end

p = [1.0, 1.5, 0.1, 0.01, 0.03, 1500.0, 0.0]
t = (0.0, 2.0)
u0 = [0.0, -0.1, 0.0]

prob = ODEProblem(single_stage, u0, t, p)
sol = solve(prob)
plot(sol)

#objective function to optimize coil gun
#let's vary: cutoff_time
function ob!(g, x)
    #set up variables
    p = [1.0, 1.5, 0.1, 0.01, 0.03, 1500.0, x[1]] #important variables
    t = (0.0, 2.0) #time
    u0 = [0.0, -0.1, 0.0] #initial conditions
    
    #simulate
    prob = ODEProblem(single_stage, u0, t, p) #set up the problem
    sol = solve(prob) #solve the problem
    solution = [u for (u) in tuples(sol)] #Extract the data into a format we can use more easily
        
    #objective (velocity when x = pos)
    pos = 0.1
    indx = find_from_sol(solution, pos) #find the index of the solution where the position of the bullet is just past the determined pos
    #make a vector of interpolated values for the position



    #changing the step size here changes the answer more than I think it should
    t = range(start = sol.t[indx-1], stop = sol.t[indx], step = 0.001) #make a range of the times from sol that bound our desired pos
    

    
    location = [sol(t)[2] for (t) in t] #make a vector of locations corresponding to t
    second_indx = find_from_loc(location, pos) #find the index of the time where the location is just greater than pos
    time = t[second_indx] #get what time it was at that index
    f = -sol(time)[3] #plug that time into sol and retreive the velocity (in effect, get the velocity of the bullet at pos)

    #constraints
    g[1] = x[1] #the voltage can't be turned off at a negative value nor one larger than 2s
    
    #return
    return f
end

function optimize_coil_gun
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
end