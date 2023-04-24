using InternalForces 


L = 600.0 
E = 29000.0
w = 0.001
I = 10.0

z = range(0.0, L, 20) 

deflection(z, L, E, I) = w * z/(24 * E * I) * (L^3 - 2 * L * z^2 + z^3)

Δ = deflection.(z, L, E, I)

M = InternalForces.moment(z, -Δ, E, I)

M_max = maximum(M)

isapprox(M_theor, w*L^2/8, rtol=0.01)

