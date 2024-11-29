using Pkg
Pkg.activate(dirname(@__FILE__))

using GLMakie

using ElectricFields
using StrongFieldApproximation

using Unitful
using UnitfulAtomic

using FFTW
import DSP: hanning

using DelimitedFiles
