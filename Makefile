CXX      = /opt/local/bin/mpic++
CC       = /usr/local/bin/gcc-4.8

CXXFLAGS = -O3 -std=c++11 -I/opt/local/include
CFLAGS   = -O3
LDFLAGS  = -L/opt/local/lib -lexoIIv2c -lnetcdf -lhdf5_hl -lhdf5  -lz -lcurl

OBJECTS		= \
	./obj/main.o \
	./obj/exodus_file.o \
	./obj/ses3d.o \
	./obj/model_file.o \
	./obj/utilities.o \
	./obj/facet.o \
	./obj/kdtree.o
	
./obj/%.o: ./src/%.cpp
	@$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "CXX $<"
	
./obj/%.o: ./src/%.c
	@$(CC) $(CFLAGS) -c -o $@ $<
	@echo "CC  $<"
	
#######
all: main

main: $(OBJECTS)
	@$(CXX) $(LDFLAGS) -o ./bin/test.exe $(OBJECTS) $(LDFLAGS)
	@echo "Linking... $<"
	
#######
clean:
	$(RM) ./obj/*.o ./bin/*.exe
