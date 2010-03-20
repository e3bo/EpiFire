bin_PROGRAMS = dynet
dynet_SOURCES = main.cpp uthash.h Network.cpp Network.h Utility.cpp Utility.h MersenneTwister.h 
LDADD = -lgsl -lgslcblas

check_PROGRAMS = debug_t1 debug_t2 debug_t3 debug_t4 debug_t5 bugsearch
debug_t1_SOURCES = $(dynet_SOURCES)
debug_t1_CPPFLAGS = -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=1
debug_t2_SOURCES = $(dynet_SOURCES)
debug_t2_CPPFLAGS = -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=2
debug_t3_SOURCES = $(dynet_SOURCES)
debug_t3_CPPFLAGS = -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=3
debug_t4_SOURCES = $(dynet_SOURCES)
debug_t4_CPPFLAGS = -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=4
debug_t5_SOURCES = $(dynet_SOURCES)
debug_t5_CPPFLAGS = -DDYNET_DEBUG=2 -DDYNET_DEBUG_TEST=5
bugsearch_SOURCES = $(dynet_SOURCES)
bugsearch_CPPFLAGS = -DDYNET_DEBUG=2 
TESTS = debug_t1 debug_t2 debug_t3 debug_t4 debug_t5 bugsearch
XFAIL_TESTS = debug_t1 debug_t2 debug_t3 debug_t4 debug_t5
