
using FFTW

# global variables
N_POINTS_P_AXIS = 50
ν = 1.0 / 1_600
TIME_STEP_LENGTH = 0.02
N_TIME_STEPS = 1_700
PLOT_EVERY = 3

function cross_product!(
    res_x::Array{T,3},
    res_y::Array{T,3},
    res_z::Array{T,3},
    a_x::Array{T,3},
    a_y::Array{T,3},
    a_z::Array{T,3},
    b_x::Array{T,3},
    b_y::Array{T,3},
    b_z::Array{T,3},
) where {T<:Union{Float64,Complex{Float64}}}
    res_x .= (a_y .* b_z .- a_z .* b_y)
    res_y .= (a_z .* b_x .- a_x .* b_z)
    res_z .= (a_x .* b_y .- a_y .* b_x)
end

function step(x)
    # (1) Compute the Curl in Fourier Domain
    cross_product!(
        curl_x_fft,
        curl_y_fft,
        curl_z_fft,
        im .* wavenumbers_x,
        im .* wavenumbers_y,
        im .* wavenumbers_z,
        velocity_x_fft,
        velocity_y_fft,
        velocity_z_fft,
    )

    # (2) Transform curl to spatial domain
    curl_x .= real(fft_operator \ curl_x_fft)
    curl_y .= real(fft_operator \ curl_y_fft)
    curl_z .= real(fft_operator \ curl_z_fft)

    # (3) Compute "Convection" in spatial domain
    cross_product!(
        convection_x, convection_y, convection_z,
        velocity_x, velocity_y, velocity_z,
        curl_x, curl_y, curl_z,
    )

    # (4) Transform "Convection" to Fourier Domain
    convection_x_fft .= fft_operator * convection_x
    convection_y_fft .= fft_operator * convection_y
    convection_z_fft .= fft_operator * convection_z

    # (5) Dealiasing on the higher wavenumbers
    convection_x_fft .*= dealias
    convection_y_fft .*= dealias
    convection_z_fft .*= dealias

    # (6) Compute the Pseudo-Pressure by a Divergence in Fourier Domain
    pressure_fft = (
        normalized_wavenumbers_x .* convection_x_fft
        + normalized_wavenumbers_y .* convection_y_fft
        + normalized_wavenumbers_z .* convection_z_fft
    )

    # (7) Assemble the rhs to the ODE system in Fourier Domain
    rhs_x_fft = (
        convection_x_fft
        -
        ν * wavenumbers_norm .^ 2 .* velocity_x_fft
        -
        normalized_wavenumbers_x .* pressure_fft
    )
    rhs_y_fft = (
        convection_y_fft
        -
        ν * wavenumbers_norm .^ 2 .* velocity_y_fft
        -
        normalized_wavenumbers_y .* pressure_fft
    )
    rhs_z_fft = (
        convection_z_fft
        -
        ν * wavenumbers_norm .^ 2 .* velocity_z_fft
        -
        normalized_wavenumbers_z .* pressure_fft
    )

    # (8) Euler Step Update
    velocity_x_fft .+= rhs_x_fft * TIME_STEP_LENGTH
    velocity_y_fft .+= rhs_y_fft * TIME_STEP_LENGTH
    velocity_z_fft .+= rhs_z_fft * TIME_STEP_LENGTH

    # (9) Transform the velocities back to spatial domain
    velocity_x .= real(fft_operator \ velocity_x_fft)
    velocity_y .= real(fft_operator \ velocity_y_fft)
    velocity_z .= real(fft_operator \ velocity_z_fft)
end
