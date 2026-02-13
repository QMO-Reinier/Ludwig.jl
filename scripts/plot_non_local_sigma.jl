using HDF5
using PlotlyJS
using ProgressMeter
using Statistics

# Open the datafile
datafile = joinpath(@__DIR__, "..", "data/conductivity", "sigma_14K_5n_e_20_n_theta_no_imp.h5")
data = h5open(datafile, "r")

# Extract the data (all have units of eV)
q_array_SI     = read(data["data/q_array"])
sigma_q_100_xx = read(data["data/sigma_q_100_xx"])
sigma_q_110_xx = read(data["data/sigma_q_110_xx"])
sigma_q_100_xy = read(data["data/sigma_q_100_xy"])
sigma_q_110_xy = read(data["data/sigma_q_110_xy"])
sigma_q_100_yx = read(data["data/sigma_q_100_yx"])
sigma_q_110_yx = read(data["data/sigma_q_110_yx"])
sigma_q_100_yy = read(data["data/sigma_q_100_yy"])
sigma_q_110_yy = read(data["data/sigma_q_110_yy"])

plt = plot(
    [
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_100_xx)/real(sigma_q_100_xx[1]),
            mode = "lines",
            name = "σ_xx q || 100",
            line = attr(color="green")
        ),
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_110_xx)/real(sigma_q_110_xx[1]),
            mode = "lines",
            name = "σ_xx q || 110",
            line = attr(color="blue")
        ),
    ],
    Layout(
        xaxis = attr(title="q (um^-1)"),
        yaxis = attr(title="σ(q)/σ(0)"),
        font = attr(size=20),
        #legend=attr(title="Bands"),
    )
)

display(plt)

plt = plot(
    [
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_100_yy)/real(sigma_q_100_yy[1]),
            mode = "lines",
            name = "σ_yy q || 100",
            line = attr(color="green")
        ),
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_110_yy)/real(sigma_q_110_yy[1]),
            mode = "lines",
            name = "σ_yy q || 110",
            line = attr(color="blue")
        ),
    ],
    Layout(
        xaxis = attr(title="q (um^-1)"),
        yaxis = attr(title="σ(q)/σ(0)"),
        font = attr(size=20),
        #legend=attr(title="Bands"),
    )
)

display(plt)

plt = plot(
    [
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_100_xx)/real(sigma_q_100_xx[1]),
            mode = "lines",
            name = "σ_xx q || 100",
            line = attr(color="green")
        ),
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_100_yy)/real(sigma_q_100_yy[1]),
            mode = "lines",
            name = "σ_yy q || 100",
            line = attr(color="blue")
        ),
    ],
    Layout(
        xaxis = attr(title="q (um^-1)"),
        yaxis = attr(title="σ(q)/σ(0)"),
        font = attr(size=20),
        #legend=attr(title="Bands"),
    )
)

display(plt)

plt = plot(
    [
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_110_xx)/real(sigma_q_110_xx[1]),
            mode = "lines",
            name = "σ_xx q || 110",
            line = attr(color="green")
        ),
        scatter(
            x = q_array_SI*1e-6,
            y = real(sigma_q_110_yy)/real(sigma_q_110_yy[1]),
            mode = "lines",
            name = "σ_yy q || 110",
            line = attr(color="blue")
        ),
    ],
    Layout(
        xaxis = attr(title="q (um^-1)"),
        yaxis = attr(title="σ(q)/σ(0)"),
        font = attr(size=20),
        #legend=attr(title="Bands"),
    )
)

display(plt)