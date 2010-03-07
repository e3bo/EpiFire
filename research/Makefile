
all: tmp_research dynet

tmp_research: tmp_research.cpp
	g++ -I ../src ../src/*.o  tmp_research.cpp -o tmp_research

dynet: dynet.cpp
	g++ -I ../src ../src/*.o  -g -Wall dynet.cpp -o dynet




