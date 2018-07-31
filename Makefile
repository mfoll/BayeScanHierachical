# BayeScan makefile
# Compile using docker: docker run --rm -v "$PWD":/usr/src/myapp -w /usr/src/myapp gcc:4.9 make

bayescanH: start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 
	g++ -O3 -fopenmp -static -o bayescanH start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 

start.o: start.cpp errors.cpp anyoption.h global_defs.h
	g++ -O3 -fopenmp -c start.cpp errors.cpp 

beta.o: beta.cpp global_defs.h
	g++ -O3 -fopenmp -c beta.cpp 
      
dirichlet.o: dirichlet.cpp global_defs.h
	g++ -O3 -fopenmp -c dirichlet.cpp 

RJupdates.o: RJupdates.cpp global_defs.h
	g++ -O3 -fopenmp -c RJupdates.cpp 

MHupdates.o: MHupdates.cpp global_defs.h
	g++ -O3 -fopenmp -c MHupdates.cpp 

likelihood.o: likelihood.cpp global_defs.h
	g++ -O3 -fopenmp -c likelihood.cpp 

read_write.o: read_write.cpp errors.cpp global_defs.h
	g++ -O3 -fopenmp -c read_write.cpp errors.cpp 

anyoption.o: anyoption.cpp anyoption.h 
	g++ -O3 -fopenmp -c anyoption.cpp 

clean: 
	rm *.o bayescanH
