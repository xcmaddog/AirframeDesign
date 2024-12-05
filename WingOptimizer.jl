using VortexLattice
using SNOW
using LinearAlgebra

"""
get_filter_matrix(size, R, spacing)
This function builds and returns a filter matrix for the optimizer.

size is the number of chord segments
R is the radius
spacing is the distance between each chord segment

W is the filter matrix such that xPrime = W * x where x are the original 
    variables chosen by the optimizer
"""
function get_filter_matrix(size, R, spacing)
    w = zeros(Float64, size, size)
    W = zeros(Float64, size, size)

    #build w
    for i in range(start = 1, stop = size)
        for j in range(start = 1, stop = size)
            w[i, j] = max(0, R - abs((i-j)*spacing))
        end
    end
    #build W
    for i in range(start = 1, stop = size)
        for j in range(start = 1, stop = size)
            denominator = 0
            for k in range(start = 1, stop = size)
                denominator = denominator + w[i,k]
            end
            W[i, j] = w[i, j] / denominator
        end
    end
    return W
end

"""
build_grid(lengths)
This function takes a vector of chord lengths 
    and returns the points to plug into VortexLattice 
"""
#I probably won't use this function, but I already wrote it,
#so I'm not throwing it away yet
function build_grid(chord_lengths)
    #calculate the leading and trailing positions
    leading_positions = chord_lengths .* (-1/4)
    trailing_positions = chord_lengths .* (3/4)
    y_positions = range(start = 0, stop = 4, length = length(chord_lengths))
    #initialize the grid
    grid = zeros(Float64, 3, 2, length(chord_lengths))
    #fill the grid
    for i in 1:length(chord_lengths) #traverse spanwise
        #define the x and y for the leading position
        grid[1, 1, i] = leading_positions[i]
        grid[2, 1, i] = y_positions[i]
        #define the x and y for the trailing position
        grid[1, 2, i] = trailing_positions[i]
        grid[2, 2, i] = y_positions[i]
    end
    return grid
end

"""
build_wing(chord_lengths)
This function takes a vector of chord lengths, 
    builds a wing aligning the quarter chords,
    and returns that wing along with its surface area
"""
function build_wing(chord_lengths, wing_span)
    #collect the necessary variables
    xle = chord_lengths .* (-1/4) #leading edge x-coordinate of each airfoil section
    yle = collect(range(start = 0, stop = wing_span/2, length = length(chord_lengths))) #leading edge y-coordinate of each airfoil section
    zle = zeros(Float64, length(chord_lengths)) #leading edge z-coordinate of each airfoil section
    theta = (0.001*(pi/180)) * ones(Float64, length(chord_lengths)) #twist of each airfoil section
    phi = zeros(Float64, length(chord_lengths)) #dihedral angle of each airfoil section
    ns = length(chord_lengths) - 1 #number of spanwise panels
    nc = 9 #number of chordwise panels
    #build the wing
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord_lengths, theta, phi, ns, nc)
    #calculate the area and part of the mean aerodynamic chord
    reference_area = 0.0
    chord_sum = 0.0
    y_d = (4/length(chord_lengths))
    for i = 1:length(chord_lengths)-1
        y_1 = i * y_d
        y_0 = y_1 - y_d
        reference_area = reference_area + (1/2) * y_d * (chord_lengths[i] + chord_lengths[i+1])
        chord_sub_sum = (y_d) * (chord_lengths[i] + y_0) + ((chord_lengths[i+1]) / (2*y_d)) * ((y_0^2) - (y_1^2))
        chord_sum = chord_sum + chord_sub_sum
    end
    #calculate the mean aerodynamic chord
    reference_chord = (2/reference_area) * chord_sum

    return surface, reference_area, reference_chord
end

function CL_to_L(lift_coefficient, reference_area, rho, freestream_velocity)
    return (1/2) * rho * (freestream_velocity^2) * reference_area * lift_coefficient
end

function CD_to_D(drag_coefficient, reference_area, rho, freestream_velocity)
    return (1/2) * rho * (freestream_velocity^2) * reference_area * drag_coefficient
end

"""
simulation(wing_surface, reference_area, reference_chord)
This function takes in some minimal descriptors of 
    the wing and runs a generic simulation on it.
    It returns the coefficients of lift and drag.
"""
function simulation(wing_surface, reference_area, reference_chord)
    # Set up freestream parameters
    alpha = 5.0 * pi / 180 # angle of attack
    beta = 0.0 # sideslip angle
    Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
    Vinf = 1.0 # reference velocity
    fs = Freestream(Vinf, alpha, beta, Omega)
    # Initialize reference location and symmetric option
    cg = 0.0
    rref = [cg, 0.0, 0.0] #
    symmetric = true
    # reference parameters (determined by the wing)
    ref = Reference(reference_area, reference_chord, 8, rref, Vinf)
    # create vector containing the surface and perform steady state analysis
    surfaces = [wing_surface]
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF # CD = drag, CY = side forces, CL = lift
    #get the properties
    properties = get_surface_properties(system)
    return CL, CD, properties
end

"""
optimize_wing(num_of_segments; filter_radius = 8/6, wing_span = 8)
This function...

filter_radius
"""
function optimize_wing(num_of_segments, W; filter_radius = 8/6, wing_span = 8)
    #define the objective function
    function objective!(g, x)
        #build the variables
        xPrime = W * x
        surface, reference_area, reference_chord = build_wing(xPrime, wing_span)
        #run the simulation
        CL, CD, properties = simulation(surface, reference_area, reference_chord)
        #objective
        rho = 1
        freestream_velocity = 1
        D = CD_to_D(CD, reference_area, rho, freestream_velocity)
        #constraints
        g[1] = CL_to_L(CL, reference_area, rho, freestream_velocity)
        #return
        return D
    end

    x0 = 0.1 * ones(Float64, num_of_segments)  # starting point
    lx = zeros(Float64, num_of_segments)  # lower bounds on x
    ux = 3.0 * ones(Float64, num_of_segments)  # upper bounds on x
    ng = 1  # number of constraints
    lg = [1.7]  # lower bounds on g
    ug = [Inf64]  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    xopt, fopt, info = minimize(objective!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)

    return xopt
end

function optimize_and_visualize(num_of_segments; filter_radius = 8/6, wing_span = 8)
    # build the filter radius
    W = get_filter_matrix(num_of_segments, filter_radius, wing_span/num_of_segments)
    #run the optimization
    segment_lengths = optimize_wing(num_of_segments, W; filter_radius = filter_radius, wing_span = wing_span)
    #draw the optimization:
    #build the variables
    surface, reference_area, reference_chord = build_wing((W^-1) *segment_lengths, 8)
    #run the simulation
    _, _, properties = simulation(surface, reference_area, reference_chord)
    write_vtk("optimizedWing", [surface], properties, symmetric = [true])
end

function somePlot(chord_lengths, wing_span)
    W = get_filter_matrix(length(chord_lengths), 8/6, 8)
    xPrime = W * chord_lengths
    surface, reference_area, reference_chord = build_wing(xPrime, wing_span)
    #run the simulation
    CL, CD, properties = simulation(surface, reference_area, reference_chord)
    write_vtk("initialWing", [surface], properties, symmetric = [true])
    #objective
    rho = 1
    freestream_velocity = 1
    D = CD_to_D(CD, reference_area, rho, freestream_velocity)
    #constraints
    L = CL_to_L(CL, reference_area, rho, freestream_velocity)
    println(CD)
    println(CL)
    println(D)
    println(L)
end

#optimize_and_visualize(10)
#x0 = 0.1 * ones(Float64, 10)  # starting point
#x0 = [1 2 3 4 5 6 7 8 9 10]'
x0 = [6.893416746978425e-9, 3.000000023674123, 3.000000017746777, 1.4998569274953353, 0.21886481079452977, 1.2022333013459494, 0.4399460841624301, 0.8209073037334909, 1.2238802596959992, 0.6405113570054716]
somePlot(x0, 8)