include Makefile.local


J123SOURCES =  bond.C bond.h solver.h J123rules.h bravaislattices.h makefile Makefile.local globalheader.h RunParameter.h



%.run  : %.exec
	./$<

%.exec : %.o 
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)


j123.o  : $(J123SOURCES)
	$(CCC) $(CCFLAGS) -DSIMPLECUBIC -DRECIPROCAL_LATTICE -c -o $@ $<        


spotless:       
	make clean
	rm -f *.ps

clean   :
	rm -f *.dvi
	rm -f *.aux
	rm -f *.log
	rm -f *~
	rm -f core
	rm -f *.o	
	rm -f *.exec
	rm -f *.d

