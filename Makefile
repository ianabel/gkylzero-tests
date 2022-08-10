
all: TaylorSedovSolution

TaylorSedovSolution: TaylorSedovSolution.c TaylorSedov.h
	$(CXX) -fpermissive -o TaylorSedovSolution -std=c17 -g -O0 TaylorSedovSolution.c -lm -lgsl
