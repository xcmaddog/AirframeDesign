using VortexLattice
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
function build_wing(chord_lengths)
    #collect the necessary variables
    xle = chord_lengths .* (-1/4) #leading edge x-coordinate of each airfoil section
    yle = collect(range(start = 0, stop = 4, length = length(chord_lengths))) #leading edge y-coordinate of each airfoil section
    zle = zeros(Float64, length(chord_lengths)) #leading edge z-coordinate of each airfoil section
    theta = (5*(pi/180)) * ones(Float64, length(chord_lengths)) #twist of each airfoil section
    phi = zeros(Float64, length(chord_lengths)) #dihedral angle of each airfoil section
    ns = length(chord_lengths) - 1 #number of spanwise panels
    sc = 10 #number of chordwise panels
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
    alpha = 0.0 * pi / 180 # angle of attack
    beta = 0.0 # sideslip angle
    Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
    Vinf = 1.0 # reference velocity
    fs = Freestream(Vinf, alpha, beta, Omega)
    # Initialize reference location and symmetric option
    cg = 0.0
    rref = SVector(cg, 0.0, 0.0) # using StaticArrays for fixed-size vector
    symmetric = true
    # reference parameters (determined by the wing)
    ref = Reference(reference_area, reference_chord, 8, rref, Vinf)
    # create vector containing the surface and perform steady state analysis
    surfaces = [wing_surface]
    system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)
    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF # CD = drag, CY = side forces, CL = lift
    return CL, CD
end

