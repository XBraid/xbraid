using ForwardDiff, Interpolations

const n_x = 8
x = range(0, 2*π*n/(n+1), n_x)
coords = zeros(3*n_x^3)
u = zeros(3*n_x^3)

function get_views(u)
    @views begin
        u_x = u[1:n_x^3]
        u_y = u[n_x^3+1:2*n_x^3]
        u_z = u[2*n_x^3+1:end]
    end
    return [u_x, u_y, u_z]
end

function advect_semi_lagrangian(u, coords)

end