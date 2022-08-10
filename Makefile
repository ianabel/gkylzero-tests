
all: TaylorSedovSolution

TaylorSedovSolution: TaylorSedovSolution.c TaylorSedov.h
	$(CC) -Wall -o TaylorSedovSolution -std=c17 -g -O0 TaylorSedovSolution.c -lm -lgsl

clean:
	rm TaylorSedovSolution;

.PHONY: clean
