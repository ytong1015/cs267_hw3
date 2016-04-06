rm *.out
module load bupc
#upcrun -n 1 -shared-heap=3G ./pgen test
upcrun -n 4 -shared-heap=1G ./pgen test
