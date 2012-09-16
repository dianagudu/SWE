
CPPFLAGS:=`pkg-config netcdf --cflags` -DCONFIG_COMPILE_WITHOUT_SIERPI=1 -DCONFIG_ENABLE_NETCDF=1 -DSIMULATION_TSUNAMI_ZERO_THRESHOLD=0.0001 -DCONFIG_DEFAULT_FLOATING_POINT_TYPE=float -I./

LDFLAGS:=`pkg-config netcdf --libs` -lnetcdf_c++


all:
	g++ -c ../../datasets_common/CSimpleNetCDF.cpp -o CSimpleNetCDF.o $(CPPFLAGS)
	g++ -c main.cpp -o main.o $(CPPFLAGS)
	g++ main.o CSimpleNetCDF.o -o main $(LDFLAGS)


clean:
	rm -f *.o
	rm -f main
