output = "q1d.out"
compile:
	gfortran types.f constants.f inputs.f setup.f utils.f common_calcs.f flx_schemes.f bc.f timestepping.f main.f -o $(output)

clean:
	rm -fv *.out *.mod *.MOD
