
all: dynet debug_dynet debug_dynet-t1 debug_dynet-t2 debug_dynet-t3 debug_dynet-t4 debug_dynet-t5

debug_dynet-t1: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=1 dynet.cpp -o debug_dynet-t1


debug_dynet-t2: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=2 dynet.cpp -o debug_dynet-t2

debug_dynet-t3: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=3 dynet.cpp -o debug_dynet-t3

debug_dynet-t4: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=4 dynet.cpp -o debug_dynet-t4

debug_dynet-t5: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=5 dynet.cpp -o debug_dynet-t5

debug_dynet: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas -DDYNET_DEBUG=2 dynet.cpp -o debug_dynet

dynet: dynet.cpp
	g++ -I ../src ../src/*.o  -g -lgsl -lgslcblas dynet.cpp -o dynet




