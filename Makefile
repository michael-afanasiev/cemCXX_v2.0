CXXFLAGS = -O3 -std=c++11 -I/Users/michaelafanasiev/Development/include
CFLAGS   = -O3
LDFLAGS  = -L/Users/michaelafanasiev/Development/lib -lexodus -lnetcdf

OBJECTS		= \
	./obj/main.o \
	./obj/exodus_file.o \
	./obj/ses3d.o \
	./obj/model_file.o
	
./obj/%.o: ./src/%.cpp
	@$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo "CXX $<"
	
#######
all: main

main: $(OBJECTS)
	@$(CXX) $(LDFLAGS) -o ./bin/test $(OBJECTS) $(LDFLAGS)
	@echo "Linking... $<"
	
#######
clean:
	$(RM) ./obj/* ./bin/*
