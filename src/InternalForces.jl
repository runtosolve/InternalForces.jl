module InternalForces

using FornbergFiniteDiff

export moment, shear, torsion, bimoment, calculateDerivativeOperators




function calculateDerivativeOperators(z)

    num_nodes = length(z)

    #first derivative
    nth_derivative = 1

    Az = zeros(Float64, (num_nodes, num_nodes))
    for i=3:num_nodes-2
    
       x = z[i-2:i+2]
       x0 = x[3]
       stencil = FornbergFiniteDiff.calculate_weights(nth_derivative, x0, x)
       Az[i, i-2:i+2] = stencil
 
    end


    #update boundary condition with single sided stencils
    order=1

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Az[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Az[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Az[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Az[end,end-4:end]=coeffs

   

    #second derivative
    Azz = zeros(Float64, (num_nodes, num_nodes))
    nth_derivative = 2
    
    for i=3:num_nodes-2
    
       x = z[i-2:i+2]
       x0 = x[3]
       stencil = FornbergFiniteDiff.calculate_weights(nth_derivative, x0, x)
       Azz[i, i-2:i+2] = stencil
 
    end


    order=2   #update stencil to consider boundary conditions without ghost nodes

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Azz[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Azz[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Azz[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Azz[end,end-4:end]=coeffs

    #third derivative

    Azzz = zeros(Float64, (num_nodes, num_nodes))
    nth_derivative = 3
    
    for i=3:num_nodes-2
    
       x = z[i-2:i+2]
       x0 = x[3]
       stencil = FornbergFiniteDiff.calculate_weights(nth_derivative, x0, x)
       Azzz[i, i-2:i+2] = stencil
 
    end


    order=3

    x=z[1:5]
    x0=x[1]
    coeffs=calculate_weights(order, x0, x)
    Azzz[1,1:5]=coeffs

    x=z[1:5]
    x0=x[2]
    coeffs=calculate_weights(order, x0, x)
    Azzz[2,1:5]=coeffs

    x=z[end-4:end]
    x0=x[end-1]
    coeffs=calculate_weights(order, x0, x)
    Azzz[end-1,end-4:end]=coeffs

    x=z[end-4:end]
    x0=x[end]
    coeffs=calculate_weights(order, x0, x)
    Azzz[end,end-4:end]=coeffs

    # #place singled sided stencils at discontinuities
    # Azzz = operatorjumps(z, dm, Azzz, order)

    return Az, Azz, Azzz

end


function moment(z, Δ, E, I)

    Az, Azz, Azzz = calculateDerivativeOperators(z)

    M = E .* I .* Azz * Δ

    return M

end

function shear(z, Δ, E, I)

    Az, Azz, Azzz = calculateDerivativeOperators(z)

    V = E .* I .* Azzz * Δ

    return V

end

function torsion(z, ϕ, E, G, J, Cw)

    Az, Azz, Azzz = calculateDerivativeOperators(z)

    #T=GJ*dϕ/dz - ECw*dϕ3/dz3
    T=G.*J.*Az*ϕ .- E.*Cw.*Azzz*ϕ

    return T

end

function bimoment(z, ϕ, E, Cw)

    Az, Azz, Azzz = calculateDerivativeOperators(z)

    #B=ECwϕ''
    B=E.*Cw.*Azz*ϕ

    return B

end

end # module
