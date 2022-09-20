#!/usr/bin/env julia
# run example to compile call to Eirene
include("xyz2PH.jl")
# remove output of example
rm("xyz2PH_test/"; recursive=true)
# write compiled image
using PackageCompiler: create_sysimage
create_sysimage(["Eirene"]; sysimage_path="xyz2PH.so")

