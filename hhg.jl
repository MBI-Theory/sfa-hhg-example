include("includes.jl")

# Find index of element e in vector v
ind(v,e) = argmin(i -> abs(v[i] - e), eachindex(v))

@field(F) do
    I₀ = 3e14u"W/cm^2"
    λ = 800.0u"nm"
    env = :tophat
    ramp = 0.0 # Length of ramps, in cycles
    flat = 5.0 # Length of flat region, in cycles
end

# # Alternative, Gaussian field (truncated smoothly from σoff to σmax,
# # given in standard deviations).
# @field(F) do
#     I₀ = 3e14u"W/cm^2"
#     λ = 800.0u"nm"
#     τ = 8u"fs"
#     env = :trunc_gauss
#     σoff = 3.0
#     σmax = 4.0
# end

# 1 electronvolt in Hartrees
eV = auconvert(u"eV", 1.0)

# Hydrogen ionization potential
Iₚ = 0.5

# Ponderomotive potential in electronvolts
Uₚ = ponderomotive_potential(F) |> u"eV"
cutoff_eV = 3.17Uₚ + Iₚ*eV
T = austrip(period(F))

# Number of time steps per cycle of the driving field
ndt = 300

t = timeaxis(F, ndt)
tplot = t*ustrip(auconvert(u"fs", 1.0))

Δt = step(t)

freq = fftshift(fftfreq(length(t), 1/Δt))
ω = 2π*freq
ω_eV = ω*ustrip(eV)

# Make a selection of the energy axis for plotting
sel = ind(ω_eV, 0):ind(ω_eV, 100)

Fv = field_amplitude(F, t)
Av = vector_potential(F, t)

# Compute the time-dependent dipole moment using the SFA, limiting the
# excursion integral to 0.65T and 2.0T, repectively, where T is the
# period time.
#
# ds and dsl are vectors with the same number of elements as our time
# grid t.
ds = induced_dipole(Iₚ, F, ndt, memory=floor(Int, 0.65T/Δt))
dsl = induced_dipole(Iₚ, F, ndt, memory=floor(Int, 2.0T/Δt))

# If you prefer to plot and analyze the data in another
# program/language, you can write the data to a text file like this:
writedlm("hhg.txt", [t ds dsl])

# Compute the Fourier transform of the induced dipole moment.
#
# It is important that you use a window function to smoothly turn off
# the signal at the beginning and the end of the time interval. Why?
# What happens if you set wind = 1 ?
wind = hanning(length(t))
Ds = fftshift(fft(ds .* wind))
Dsl = fftshift(fft(dsl .* wind))

DGs = zeros(ComplexF64, length(t), length(ω))
DGsl = zeros(ComplexF64, length(t), length(ω))

# TODO: Compute your Gabor transform here and store in the arrays DGs
# and DGsl!
for (i,τ) in enumerate(t)
    DGs[i,:] .= 0 # Replace 0 by your implementation of the Gabor transform for the time τ!
    DGsl[i,:] .= 0 # Replace 0 by your implementation of the Gabor transform for the time τ!
end

fig = Figure(size=(1200,1000))
axF = Axis(fig[1,1], xticksvisible=false, xticklabelsvisible=false, ylabel=L"$F(t)$ [au]")
axA = Axis(fig[2,1], xticksvisible=false, xticklabelsvisible=false, ylabel=L"$A(t)$ [au]")
axd = Axis(fig[3,1], xlabel=L"$t$ [fs]")

linkxaxes!(axF, axA, axd)

lines!(axF, tplot, Fv)
lines!(axA, tplot, Av)
lines!(axd, tplot, ds, label="Short trajectory only")
lines!(axd, tplot, dsl, label="Short and long trajectories")
axislegend(axd)


axD = Axis(fig[3,2], yaxisposition=:right, xlabel=L"$\omega$ [eV]", ylabel=L"$|d(\omega)|^2$ [arb.u.]",
           yscale=log10, limits=(nothing, (1e-14,nothing)))

lines!(axD, ω_eV[sel], abs2.(Ds[sel]))
lines!(axD, ω_eV[sel], abs2.(Dsl[sel]))
vlines!(axD, ustrip(cutoff_eV), linestyle=:dash)

axGs = Axis(fig[1,2],
           xaxisposition=:top,
           yaxisposition=:right,
           xlabel=L"$t$ [fs]", ylabel=L"|d(\omega)\ast G(\omega,t)|^2")
axGsl = Axis(fig[2,2],
             xticksvisible=false, xticklabelsvisible=false,
             yaxisposition=:right,
             ylabel=L"|d(\omega)\ast G(\omega,t)|^2")

linkxaxes!(axF, axA, axd, axGs, axGsl)

vmax = 1e-4
contrast = 1e-7
vmin = contrast*vmax
heatmap!(axGs, tplot, ω_eV[sel], abs2.(DGs[:,sel]),
         colorrange=(vmin,vmax), colormap=:inferno)
heatmap!(axGsl, tplot, ω_eV[sel], abs2.(DGsl[:,sel]),
         colorrange=(vmin,vmax), colormap=:inferno)

fig
