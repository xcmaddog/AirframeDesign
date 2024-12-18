using VortexLattice
using SNOW
using LinearAlgebra
using Plots


"""
build_wing(chord_lengths, wing_span; angle_of_attack)
This function takes a vector of chord lengths, a wing_span, and optionally an angle of attack. 
    Then it builds a wing aligning the quarter chords,
    and returns that wing along with its surface area, reference_chord, and the grid representing the wing
"""
function build_wing(chord_lengths, wing_span; angle_of_attack = 5, num_chord_panels = 9)
    #collect the necessary variables
    xle = chord_lengths .* (-1/4) #leading edge x-coordinate of each airfoil section
    yle = collect(range(start = 0, stop = wing_span/2, length = length(chord_lengths))) #leading edge y-coordinate of each airfoil section
    zle = zeros(Float64, length(chord_lengths)) #leading edge z-coordinate of each airfoil section
    theta = (angle_of_attack*(pi/180)) * ones(Float64, length(chord_lengths)) #twist of each airfoil section
    phi = zeros(Float64, length(chord_lengths)) #dihedral angle of each airfoil section
    ns = length(chord_lengths) - 1 #number of spanwise panels
    nc = num_chord_panels #number of chordwise panels
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

    return surface, reference_area, reference_chord, grid
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
    It returns the coefficients of lift, drag, the surface
    properties, and the System object.
"""
function simulation(wing_surface, reference_area, reference_chord)
    # Set up freestream parameters
    alpha = 0.0 * pi / 180 # angle of attack
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
    return CL, CD, properties, system
end

"""
optimize_wing(num_of_segments; wing_span = 8)
This function takes a number of segments, and optionally a wingspan and returns the optimization for such a wing.

num_of_segments is the number of chord lengths.
wing_span is the span of the wing in meters

xopt is a vector of the optimized chord lengths
"""
function optimize_wing(num_of_segments; lift_requirement = 1.7, tolerance = 1e-6, wing_span = 8)
    #define the objective function
    function objective!(g, x)
        surface, reference_area, reference_chord, _ = build_wing(x, wing_span)
        #run the simulation
        CL, CD, properties, _ = simulation(surface, reference_area, reference_chord)
        #objective
        rho = 1
        freestream_velocity = 1
        D = CD_to_D(CD, reference_area, rho, freestream_velocity)
        #constraints
        g[1] = CL_to_L(CL, reference_area, rho, freestream_velocity)
        differences = diff(x) * -1
        for i in range(2, num_of_segments)
            g[i] = differences[i-1]
        end
        #return
        return D
    end

    x0 = 0.1 * ones(Float64, num_of_segments)  # starting point
    lx = zeros(Float64, num_of_segments)  # lower bounds on x
    ux = Inf64 * ones(Float64, num_of_segments)  # upper bounds on x
    ng = num_of_segments  # number of constraints
    lg = [lift_requirement]  # lower bounds on g
    append!(lg, zeros(num_of_segments - 1))
    ug = Inf64 * ones(num_of_segments)  # upper bounds on g
    # ----- set some options ------
    ip_options = Dict(
        "tol" => tolerance
        )
    solver = IPOPT(ip_options)
    options = Options(;solver)

    xopt, fopt, info = minimize(objective!, x0, ng, lx, ux, lg, ug, options)

    println("xstar = ", xopt)
    println("fstar = ", fopt)
    println("info = ", info)

    return xopt
end

function plotLiftDistribution(system, grid, num_span_segments, plot_name; num_chord_segments = 9, span = 8)
    r, c = lifting_line_geometry([grid])
    cf, cm = lifting_line_coefficients(system, r, c)
    lift_coefficients = cf[1][3,:]
    span_locations = range(start = span/2, stop = span, length = num_span_segments)
    p = plot(span_locations, lift_coefficients,
        ylabel = "Lift Coefficient",
        xlabel = "Spanwise Location",
        label = "Wing Lift Distribution")
    max_lift_coefficient = maximum(lift_coefficients)
    #y = (b/a)*sqrt(a^2 - x^2)
    elliptical_lift_coefficients = (max_lift_coefficient/(span/2)) .* sqrt.((span/2)^2 .- span_locations)
    plot!(p, span_locations, elliptical_lift_coefficients,
        label = "Elliptical Lift Distribution")
    savefig(p, plot_name * ".jpg")
end

function optimize_and_visualize(num_of_segments, vtk_name; lift_requirement = 1.7, tolerance = 1e-6, wing_span = 8)
    #run the optimization
    segment_lengths = optimize_wing(num_of_segments, lift_requirement = lift_requirement, tolerance = tolerance, wing_span = wing_span)
    #draw the optimization:
    #build the variables
    surface, reference_area, reference_chord, grid = build_wing(segment_lengths, wing_span)
    #run the simulation
    _, _, properties, system = simulation(surface, reference_area, reference_chord)
    #make the vtk
    write_vtk(vtk_name, [surface], properties, symmetric = [true])
    #make a plot that compares its lift distribution to an elliptical lift distribution
    #plotLiftDistribution(system, grid, num_of_segments, vtk_name * "_lift_distribution", span = wing_span)
    plotLiftDistribution(system, grid, num_of_segments, "distribution", span = wing_span)
end

function somePlot(chord_lengths, wing_span)
    surface, reference_area, reference_chord, _ = build_wing(x, wing_span)
    #run the simulation
    CL, CD, properties, _ = simulation(surface, reference_area, reference_chord)
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

optimize_and_visualize(10, "optimizedWing", lift_requirement = 1.7, tolerance = 1e-6, wing_span = 8) #this is the "default"

#optimize_and_visualize(35, "looseTol", lift_requirement = 1.7, tolerance = 1e-5, wing_span = 8) #this is loosening the tolerance
#optimize_and_visualize(35, "tightTol", lift_requirement = 1.7, tolerance = 1e-7, wing_span = 8) #this is tightening the tolerance

#optimize_and_visualize(5, "moreLiftWing", lift_requirement = 1.7*10, tolerance = 1e-6, wing_span = 8) #this is trying to make the wing give 10x the lift