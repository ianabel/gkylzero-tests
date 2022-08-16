
all: TaylorSedovTest

G0_DIR ?= /home/ian/projects/gkylzero

CFLAGS=-O3 -g -ffast-math -fPIC -Wall
SRC_DIRS := minus zero apps kernels
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB_DIR = ${HOME}/gkylsoft/OpenBLAS/lib
LAPACK_LIB = ${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a
INCLUDES = -I$(G0_DIR) -I$(G0_DIR)/minus -I$(G0_DIR)/minus/STC/include -I$(G0_DIR)/zero -I$(G0_DIR)/apps -I$(G0_DIR)/regression -I${LAPACK_INC} -I${SUPERLU_INC}

#Don't think I need the following includes
#INCLUDES = -I$(G0_DIR) -I$(G0_DIR)/zero -I$(G0_DIR)/apps -I$(G0_DIR)/regression -I$(G0_DIR)/minus -I$(G0_DIR)/minus/STC -I$(G0_DIR)/minus/STC/docs -I$(G0_DIR)/minus/STC/docs/pics -I$(G0_DIR)/minus/STC/include -I$(G0_DIR)/minus/STC/include/stc -I$(G0_DIR)/minus/STC/include/stc/alt -I$(G0_DIR)/zero -I$(G0_DIR)/apps -I$(G0_DIR)/kernels -I$(G0_DIR)/kernels/advection -I$(G0_DIR)/kernels/fem_parproj -I$(G0_DIR)/kernels/basis -I$(G0_DIR)/kernels/lbo -I$(G0_DIR)/kernels/fem_poisson -I$(G0_DIR)/kernels/maxwell -I$(G0_DIR)/kernels/vlasov -I$(G0_DIR)/kernels/dg_diffusion

# SuperLU includes and librararies
SUPERLU_INC = ${HOME}/gkylsoft/superlu/include
ifeq ($(UNAME_S),Linux)
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib64
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib64/libsuperlu.a
else
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib/libsuperlu.a
endif

TaylorSedovTest: rt_euler_taylorsedov.c TaylorSedovSolution.c TaylorSedov.h
	$(CC) $(CFLAGS) $(INCLUDES) -o TaylorSedovTest rt_euler_taylorsedov.c TaylorSedovSolution.c ${G0_DIR}/build/libgkylzero.a -lm -lgsl ${SUPERLU_LIB} ${LAPACK_LIB} -lpthread

rt_euler_taylorsedov_2d: rt_euler_taylorsedov_2d.c TaylorSedovSolution.c TaylorSedov.h
	$(CC) $(CFLAGS) $(INCLUDES) -o rt_euler_taylorsedov_2d rt_euler_taylorsedov_2d.c TaylorSedovSolution.c ${G0_DIR}/build/libgkylzero.a -lm -lgsl ${SUPERLU_LIB} ${LAPACK_LIB} -lpthread

TaylorSedovSolution.o: TaylorSedovSolution.c TaylorSedov.h
	$(CC) -c -Og -g -fPIC -o TaylorSedovSolution.o TaylorSedovSolution.c

TaylorSedovSolution.so: TaylorSedovSolution.o
	ld -shared -o TaylorSedovSolution.so TaylorSedovSolution.o -lm -lgsl

TaylorSedovMain: TaylorSedovMain.c TaylorSedovSolution.so 
	$(CC) -Og -g -fPIC -o TaylorSedovMain TaylorSedovMain.c -L. -Wl,-rpath . -l:TaylorSedovSolution.so -lm -lgsl


clean:
	rm -f TaylorSedovTest TaylorSedovSolution.so TaylorSedovSolution.o rt_euler_taylorsedov_2d *.gkyl *.json TaylorSedovMain;

.PHONY: clean
