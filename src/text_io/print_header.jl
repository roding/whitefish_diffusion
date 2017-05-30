include("get_version.jl")

function print_header()
	version_string::String = get_version()
	println(join(("This is Whitefish diffusion module version ", version_string, ", executed on Julia version ", string(Sys.VERSION), " (", Sys.MACHINE, ")")))
	
	nothing	
end
