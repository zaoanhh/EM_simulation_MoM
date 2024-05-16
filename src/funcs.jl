module MoMForward
using ExportAll
using LinearAlgebra
using MKL
using SpecialFunctions
using FFTW
@exportAll()

### MoMForward

@doc raw"""
    get_antenna_pos(R,N)

Get the position of the `N` antennas uniformly distributed on a circle with radius `R`. The result is a 2xN matrix, where the first row is the x coordinates and the second row is the y coordinates.
"""
function get_antenna_pos(R, N::Int)
    θ = 2π / N
    x = R * cos.(θ * (0:(N - 1)))
    y = R * sin.(θ * (0:(N - 1)))
    return [x y]'
end

"""
    ConstantParameter(;
    freq,
    doi_size_x,
    doi_size_y,
    grid_number_x,
    grid_number_y,
    R_tx,
    R_rx,
    tx_num,
    rxs_num
    )

Construct a structure based on the frequency of the input, including the following parameters:

- `ϵ₀` is the permittivity of free space.
- `μ₀` is the permeability of free space.
- `c` is the speed of light in free space.
- `λ` is the wavelength.
- `f` is the frequency.
- `ω` is the angular frequency.
- `k₀` is the wave number in the vacuum.
- `η` is the intrinsic impedance of free space.
- `constant1` is a constant factor, which is μ₀ / (4 * pi).
- `constant2` is a constant factor, which is -1im / (4 * pi * ω * ϵ₀).
- `doi_size_x` is the size of the DOI in the x direction.
- `doi_size_y` is the size of the DOI in the y direction.
- `grid_number_x` is the number of grids in the x direction.
- `grid_number_y` is the number of grids in the y direction.
- `grid_num` is the total number of grids, which is `grid_number_x * grid_number_y`.
- `grid_radius` is the radius of a circle with the same area as the grid.
- `grid_xs` is the x coordinates of the grids.
- `grid_ys` is the y coordinates of the grids.
- `txs_pos` is the position of the transmitters.
- `rxs_pos` is the position of the receivers.
- `tx_num` is the number of transmitters.
- `rxs_num` is the number of receivers.
# Keyword Arguments

- `freq::Float64`: The frequency of the input.
- `doi_size_x::Float64`: The size of the DOI in the x direction.
- `doi_size_y::Float64`: The size of the DOI in the y direction
- `grid_number_x::Int64`: The number of grids in the x direction.
- `grid_number_y::Int64`: The number of grids in the y direction.
- `R_tx`: The radius of the transmitter circle.
- `R_rx`: The radius of the receiver circle.
- `tx_num`: The number of transmitters.
- `rxs_num`: The number of receivers.
"""
struct ConstantParameter
    ϵ₀::Float64
    μ₀::Float64
    c::Float64
    λ::Float64
    f::Float64
    ω::Float64
    k₀::Float64
    η::Float64
    constant1::ComplexF64
    constant2::ComplexF64
    doi_size_x::Float64
    doi_size_y::Float64
    grid_number_x::Int64
    grid_number_y::Int64
    grid_num::Int64
    grid_radius::Float64
    grid_xs::Array{Float64, 1}
    grid_ys::Array{Float64, 1}
    txs_pos::Array{Float64, 2}
    rxs_pos::Array{Float64, 2}
    tx_num::Int64
    rxs_num::Int64
end
function ConstantParameter(;
        freq,
        doi_size_x,
        doi_size_y,
        grid_number_x,
        grid_number_y,
        R_tx,
        R_rx,
        tx_num,
        rxs_num
)
    f = freq
    ϵ₀ = 8.854e-12
    μ₀ = 1.257e-6
    c = 1.0 / sqrt(ϵ₀ * μ₀)
    λ = c / f
    ω = 2.0 * pi * f
    k₀ = ω / c
    η = sqrt(μ₀ / ϵ₀)
    constant1 = ComplexF64(μ₀ / (4 * pi))
    constant2 = ComplexF64(-1im / (4 * pi * ω * ϵ₀))
    @assert grid_number_x % 2 == 0&&grid_number_y % 2 == 0 "grid_number_x and grid_number_y should be even"
    @assert sqrt(doi_size_x^2 + doi_size_y^2)<2 * R_tx "The DOI should be smaller than the transmitter circle"
    @assert sqrt(doi_size_x^2 + doi_size_y^2)<2 * R_rx "The DOI should be smaller than the receiver circle"
    txs_pos = get_antenna_pos(R_tx, tx_num)
    rxs_pos = get_antenna_pos(R_rx, rxs_num)
    tx_num = size(txs_pos, 2)
    rxs_num = size(rxs_pos, 2)
    grid_size_x = doi_size_x / grid_number_x
    grid_size_y = doi_size_y / grid_number_y
    radius = sqrt(grid_size_x * grid_size_y / pi)
    grid_xs = ((-(grid_number_x - 1) / 2):1:((grid_number_x - 1) / 2)) .* grid_size_x
    grid_ys = (((grid_number_y - 1) / 2):-1:(-(grid_number_y - 1) / 2)) .* grid_size_y
    return ConstantParameter(
        ϵ₀, μ₀, c, λ, f, ω, k₀, η, constant1, constant2, doi_size_x, doi_size_y,
        grid_number_x, grid_number_y, grid_number_x * grid_number_y,
        radius, grid_xs, grid_ys, txs_pos, rxs_pos, tx_num, rxs_num)
end
function ConstantParameter(;
        freq, doi_size_x, doi_size_y, grid_number_x, grid_number_y, txs_pos, rxs_pos)
    f = freq
    ϵ₀ = 8.854e-12
    μ₀ = 1.257e-6
    c = 1.0 / sqrt(ϵ₀ * μ₀)
    λ = c / f
    ω = 2.0 * pi * f
    k₀ = ω / c
    η = sqrt(μ₀ / ϵ₀)
    constant1 = ComplexF64(μ₀ / (4 * pi))
    constant2 = ComplexF64(-1im / (4 * pi * ω * ϵ₀))
    @assert grid_number_x % 2 == 0&&grid_number_y % 2 == 0 "grid_number_x and grid_number_y should be even"
    txs_pos = txs_pos
    rxs_pos = rxs_pos
    tx_num = size(txs_pos, 2)
    rxs_num = size(rxs_pos, 2)
    grid_size_x = doi_size_x / grid_number_x
    grid_size_y = doi_size_y / grid_number_y
    radius = sqrt(grid_size_x * grid_size_y / pi)
    grid_xs = ((-(grid_number_x - 1) / 2):1:((grid_number_x - 1) / 2)) .* grid_size_x
    grid_ys = (((grid_number_y - 1) / 2):-1:(-(grid_number_y - 1) / 2)) .* grid_size_y
    return ConstantParameter(
        ϵ₀, μ₀, c, λ, f, ω, k₀, η, constant1, constant2, doi_size_x, doi_size_y,
        grid_number_x, grid_number_y, grid_number_x * grid_number_y,
        radius, grid_xs, grid_ys, txs_pos, rxs_pos, tx_num, rxs_num)
end
@doc raw"""
    find_grids_with_object(config, scatter; bacground_permittivity = 1.0)

We find the position in the matrix and the of coordinates the scatter grid, whose complex permittivity is not equal to that of the background.

# Arguments
- `config::ConstantParameter`: The structure of the configuration.
- `scatter::Array`: The permittivity of the object and background, which is a 2D array of booleans.

# Keyword Arguments
- `bacground_permittivity::Float64`: The permittivity of the background, which is 1.0 by default.

# Returns
- `scatter_positions_in_matrix::Array`: The position in the matrix of the scatter grid.
- `scatter_centroids::Array`: The coordinates of the scatter grid
"""
function find_grids_with_object(config::ConstantParameter,
        scatter;
        bacground_permittivity = 1.0)
    scatter_positions_in_matrix = findall(scatter .!= bacground_permittivity)
    scatter_centroid_xs = [config.grid_xs[scatter_positions_in_matrix[i][2]]
                           for
                           i in axes(scatter_positions_in_matrix, 1)]
    scatter_centroid_ys = [config.grid_ys[scatter_positions_in_matrix[i][1]]
                           for
                           i in axes(scatter_positions_in_matrix, 1)]
    object_grid_indices = findall(vec(scatter) .!= bacground_permittivity)
    return scatter_centroid_xs, scatter_centroid_ys, object_grid_indices
end

function get_Z_for_one_grid(config::ConstantParameter,
        grid_now_index,
        grid_centroid_all_x,
        grid_centroid_all_y,
        permittivity;)
    grid_centroid_all_x = grid_centroid_all_x
    grid_centroid_all_y = grid_centroid_all_y
    grid_now_x = grid_centroid_all_x[grid_now_index]
    grid_now_y = grid_centroid_all_y[grid_now_index]

    dipole_distances = sqrt.((grid_now_x .- grid_centroid_all_x) .^ 2 .+
                             (grid_now_y .- grid_centroid_all_y) .^ 2,)

    Z = (-config.η *
         π *
         (config.grid_radius / 2) *
         besselj(1, config.k₀ * config.grid_radius) .*
         hankelh1.(0, config.k₀ * dipole_distances))
    Z[grid_now_index] = -config.η *
                        π *
                        (config.grid_radius / 2) *
                        hankelh1(1, config.k₀ * config.grid_radius) -
                        1im * config.η * permittivity / (config.k₀ * (permittivity - 1))

    return Z
end

function get_impedance_in_domain(config::ConstantParameter,
        scatter_grids_permittivity,
        scatter_centroid_xs,
        scatter_centroid_ys;)
    m = length(scatter_centroid_xs)

    object_field = zeros(ComplexF64, m, m)

    for i in 1:m
        object_field[i, :] = get_Z_for_one_grid(config,
            i,
            scatter_centroid_xs,
            scatter_centroid_ys,
            scatter_grids_permittivity[i])
    end
    return object_field
end

function get_Ei_in_rxs(config::ConstantParameter)
    receiver_x = config.rxs_pos[1, :]
    receiver_y = config.rxs_pos[2, :]
    transmitter_x = config.txs_pos[1, :]
    transmitter_y = config.txs_pos[2, :]

    xtd = transmitter_x' #横着是tx
    xrd = receiver_x
    ytd = transmitter_y'
    yrd = receiver_y
    dist = sqrt.((xtd .- xrd) .^ 2 + (ytd .- yrd) .^ 2)
    # Ei_in_rxs = (1im / 4) .* hankelh1.(0, config.k₀ .* dist)
    Ei_in_rxs = exp.(-1im .* config.k₀ .* dist) ./ dist
    return Ei_in_rxs
end

function get_Ei_in_domain(config::ConstantParameter)
    transmitter_x = config.txs_pos[1, :]
    transmitter_y = config.txs_pos[2, :]
    grid_xs = config.grid_xs
    grid_ys = config.grid_ys
    dist = zeros(config.grid_num, config.tx_num)
    for tx in 1:(config.tx_num)
        temp_dist = vec(sqrt.((transmitter_x[tx] .- transpose(grid_xs)) .^ 2 .+
                              (transmitter_y[tx] .- grid_ys) .^ 2))
        dist[:, tx] = temp_dist
    end
    # Ei_in_domain = (1im / 4) .* hankelh1.(0, config.k₀ .* dist)
    Ei_in_domain = exp.(-1im .* config.k₀ .* dist) ./ dist
    return Ei_in_domain
end

function get_induced_current(config::ConstantParameter,
        impedance_in_domain,
        Ei_in_scatterer,
        object_grid_indices)
    J1 = ((impedance_in_domain) \ (-Ei_in_scatterer))
    current = zeros(ComplexF64, config.grid_num, config.tx_num)
    # for i in axes(object_grid_indices, 1)
    #     current[object_grid_indices[i], :] = J1[i, :]
    # end
    current[object_grid_indices, :] .= J1
    return current
end

function get_scattered_field(config::ConstantParameter,
        current,
        GS;
)
    scattered_field = (GS) * (current)
    return scattered_field
end

function get_GS(config::ConstantParameter)
    rx_xs = config.rxs_pos[1, :]
    rx_ys = config.rxs_pos[2, :]
    grid_xs = config.grid_xs
    grid_ys = config.grid_ys
    dist = zeros(config.grid_num, config.rxs_num)
    for rx in 1:(config.rxs_num)
        temp_dist = vec(sqrt.((rx_xs[rx] .- transpose(grid_xs)) .^ 2 .+
                              (rx_ys[rx] .- grid_ys) .^ 2))
        dist[:, rx] = temp_dist
    end

    GS = -config.η *
         π *
         (config.grid_radius / 2) *
         besselj(1, config.k₀ * config.grid_radius) .*
         hankelh1.(0, config.k₀ .* transpose(dist))
    return GS
end
@doc raw"""
    get_GD_Z(config::ConstantParameter)
"""
function get_GD_Z(config::ConstantParameter)
    M = config.grid_number_x
    N = config.grid_number_y
    grid_size_x = config.doi_size_x / config.grid_number_x
    grid_size_y = config.doi_size_y / config.grid_number_y
    x_diff = ((1 - M):1:(M - 1)) .* grid_size_x
    y_diff = ((1 - N):1:(N - 1)) .* grid_size_y

    R = sqrt.(transpose(x_diff .^ 2) .+ y_diff .^ 2)
    ZZ = -config.η *
         π *
         (config.grid_radius / 2) *
         besselj(1, config.k₀ * config.grid_radius) .* hankelh1.(0, config.k₀ .* R)
    ZZ[N, M] = -config.η *
               pi *
               (config.grid_radius / 2) *
               hankelh1(1, config.k₀ * config.grid_radius) - 1im * config.η / (config.k₀)

    Z = zeros(ComplexF64, 2 * N - 1, 2 * M - 1)
    Z[1:N, 1:M] = ZZ[N:(2 * N - 1), M:(2 * M - 1)]
    Z[(N + 1):(2 * N - 1), (M + 1):(2 * M - 1)] = ZZ[1:(N - 1), 1:(M - 1)]
    Z[1:N, (M + 1):(2 * M - 1)] = ZZ[N:(2 * N - 1), 1:(M - 1)]
    Z[(N + 1):(2 * N - 1), 1:M] = ZZ[1:(N - 1), M:(2 * M - 1)]

    return Z
end

@doc raw"""
    GD(J, Z)

We calculate the Green's function of the domain using the inverse Fourier transform. $G_D J$ is `GD(J, Z)`, where `Z` is given by [`get_GD_Z`](@ref).
"""
function GD(J, Z)
    Ni = size(J, 2)
    N, M = Int.((size(Z) .+ 1) ./ 2)
    J_now = zeros(eltype(J), 2 * N - 1, 2 * M - 1, Ni)
    J_now[1:N, 1:M, :] = reshape(J, N, M, Ni)
    opa = zeros(eltype(J), N * M, Ni)
    Z_fft = fft(Z)
    for tx in 1:Ni
        temp_opa = ifft(Z_fft .* fft(J_now[:, :, tx]))
        opa[:, tx] = vec(temp_opa[1:N, 1:M])
    end
    return opa
end
function get_shared_variables(parameters)
    Ei_in_rxs = get_Ei_in_rxs(parameters)
    Ei_in_domain = get_Ei_in_domain(parameters)
    GS = get_GS(parameters)
    GDZ = get_GD_Z(parameters)

    return Dict("Ei_in_rxs" => Ei_in_rxs,
        "Ei_in_domain" => Ei_in_domain,
        "GS" => GS,
        "GDZ" => GDZ)
end
end

module Simu
using ExportAll, LinearAlgebra, MKL, Random, Distributions
using ..MoMForward: ConstantParameter
import ..MoMForward as mf
using DrWatson
# function awgn(X, SNR;)
#     #Assumes X to be a matrix and SNR a signal-to-noise ratio specified in decibel (dB)
#     #Implented by author, inspired by https://www.gaussianwaves.com/2015/06/how-to-generate-awgn-noise-in-matlaboctave-without-using-in-built-awgn-function/
#     N = length(X) #Number of elements in X
#     signalPower = sum(X[:] .^ 2) / N
#     linearSNR = 10^(SNR / 10)
#     a, b = size(X)
#     noiseMat = (randn((a, b))) .* √(signalPower / linearSNR) #Random gaussian noise, scaled according to the signal power of the entire matrix (!) and specified SNR

#     return X + noiseMat
# end
# function awgn(signal::AbstractArray{T}, snr_db::Real) where T<:Number
#     snr_linear = 10^(snr_db / 10)
#     signal_power = mean(abs2.(signal))
#     noise_variance = signal_power / snr_linear
#     noise = sqrt(noise_variance) * randn(T, size(signal)...)
#     return signal + noise
# end
function get_Es(parameters::ConstantParameter, shared_variables::Dict, scatterer::Matrix)
    Ei_in_rxs = shared_variables["Ei_in_rxs"]
    Ei_in_domain = shared_variables["Ei_in_domain"]
    GS = shared_variables["GS"]
    GDZ = shared_variables["GDZ"]
    scatter_centroid_xs, scatter_centroid_ys, object_grid_indices = mf.find_grids_with_object(
        parameters,
        scatterer)
    impedance_in_domain = mf.get_impedance_in_domain(parameters,
        vec(scatterer)[object_grid_indices],
        scatter_centroid_xs,
        scatter_centroid_ys)
    current = mf.get_induced_current(parameters,
        impedance_in_domain,
        Ei_in_domain[object_grid_indices, :],
        object_grid_indices;)
    Es = mf.get_scattered_field(parameters, current, GS;)#The size is [num_RX, num_TX]
    return Es
end

function get_parameters_gaussian(config::Dict)
    @unpack freq, Rtx, Rrx, N, tx_num, rxs_num, Dsize, grids = config
    parameters = mf.ConstantParameter(
        freq = freq, doi_size_x = Dsize, doi_size_y = Dsize, grid_number_x = grids,
        grid_number_y = grids, R_tx = Rtx, R_rx = Rrx, tx_num = tx_num, rxs_num = rxs_num
    )
    return parameters
end
function sim_share_variables(parameters)
    shared_variables = mf.get_shared_variables(parameters)
    return shared_variables
end

function sim_Es_gaussian(config::Dict)
    Random.seed!(101)
    @unpack N, μr, μi, σr, σi, grids = config
    parameters = get_parameters_gaussian(config)
    shared_variables, _ = produce_or_load(
        sim_share_variables, parameters, datadir("shared_variables"); verbose = false
    )
    sr = Normal(μr, σr)
    si = Normal(μi, σi)
    Es_all = Vector{Matrix}(undef, N)
    epsilon_all = Vector{Vector}(undef, N)
    for i in 1:Nend
        temp_scatterer = zeros(ComplexF64, grids, grids)
        for j in 1:grids
            for k in 1:grids
                temp_scatterer[j, k] = complex(rand(sr), rand(si))
            end
        end
        Es = Simu.get_Es(parameters, shared_variables, temp_scatterer) #The size is [num_RX, num_TX]
        Es_all[i] = Es
        epsilon_all[i] = vec(temp_scatterer)
    end
    return @strdict Es_all epsilon_all
end

@exportAll()
end
