HDF5FLAGS = `pkg-config hdf5 --cflags --libs 2>/dev/null`
LGMFLAGS = `pkg-config lgm --cflags --libs`
#GSLFLAGS = `pkg-config gsl --cflags --libs`


ptm_lanlgeomag: ptm_lanlgeomag.c
	gcc -g ptm_lanlgeomag.c $(LGMFLAGS) $(HDF5FLAGS) -o ptm_lanlgeomag -Wall

clean: 
	rm -rf ptm_lanlgeomag
	echo clean done

