CXXFLAGS 	= -O3 -std=c++11
CFLAGS   	= -O3

OBJECTS		= \
	./obj/main.o \
	./obj/exodus_file.o
	
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
